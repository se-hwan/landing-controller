%% metadata
% Description: Phase-based, end effector parametrized trajectory optimization 
%              for quadrupedal landing with single rigid body model
% Author:      Se Hwan Jeon

%% cleanup
rmpath(genpath('.')); % clear all previously added paths
clear; clc; close all;

%% flags
show_animation = true;

%% add library paths
% may need to specify OS directory
% addpath(genpath('../../utilities/casadi/casadi_windows')); 
% addpath(genpath('../../utilities/casadi/casadi_linux'));
addpath(genpath('../../utilities_general'));
addpath(genpath('utilities_eeParam'));
import casadi.*

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModelMC3D(params, model, 0);

%% problem parameters and setup

dt_dyn_val = 0.2;                                               % timestep for dynamics
dt_baseSpline = 0.2;                                            % duration of polynomials for base motion
T_val = 1.0;                                                    % total time
N_baseSpline = round(T_val/dt_baseSpline);                      % number of polynomials for base
N_timesteps = round(T_val/dt_dyn_val + 1);                      % number of steps for dynamics
phase_durations_base = dt_baseSpline*ones(1,N_baseSpline);      % durations of base splines
eps = 1e-3;                                                     % small change in time to register changes bw polynomial sets

% contact sequence specification
% num contact phases, num flight phases, order (1 to begin in contact, 0 to begin in flight)
contactPhases = [2 1 1;         
                 1 2 0;         
                 1 2 0;
                 2 1 1]; 
contactSequence = generateContactSequence(contactPhases);   % contact sequence for legs
contactState_forces = cell(model.NLEGS, 1);                 % contact state for force splines
contactState_posns = cell(model.NLEGS, 1);                  % contact state for posn splines
contactDuration_forces = contactState_forces;               % contact duration for force splines
contactDuration_posns = contactState_posns;                 % contact durations for posn splines

%% optimization - variables

opti = casadi.Opti();

% decision variables:

% CoM linear positions, 4th order polynomial parametrization, fixed duration
com_x   = opti.variable(4, N_baseSpline);
com_y   = opti.variable(4, N_baseSpline);
com_z   = opti.variable(4, N_baseSpline);

% base euler angles, 4th order polynomial parametrization, fixed duration
roll    = opti.variable(4, N_baseSpline);
pitch   = opti.variable(4, N_baseSpline);
yaw     = opti.variable(4, N_baseSpline);

% stance durations
phase_durations = cell(numel(contactSequence), 1);
for i = 1:numel(contactSequence)
    phase_durations{i} = opti.variable(1, numel(contactSequence{i}));
end

% contact states and durations of splines
for i = 1:model.NLEGS
    temp_forces = [];
    temp_posns = [];
    temp_f_dur = [];
    temp_p_dur = [];
    for j = 1:length(contactSequence{i})
        if (contactSequence{i}(j) == 0)
            % not in contact
            temp_forces = [temp_forces 0];
            temp_posns = [temp_posns 0 0];
            temp_f_dur = [temp_f_dur phase_durations{i}(j)];
            temp_p_dur = [temp_p_dur ones(1,2)*(phase_durations{i}(j)/2)];
        else
            % in contact
            temp_forces = [temp_forces ones(1, 3)];
            temp_posns = [temp_posns 1];
            temp_f_dur = [temp_f_dur ones(1,3)*(phase_durations{i}(j)/3)];
            temp_p_dur = [temp_p_dur phase_durations{i}(j)];
        end
    end
    contactState_forces{i} = temp_forces;
    contactState_posns{i} = temp_posns;
    contactDuration_forces{i} = temp_f_dur;
    contactDuration_posns{i} = temp_p_dur;
end


% for end effector positions and forces, coefficients will be ordered as 
% [x0, x0_dot, x1, x1_dot]' for their coefficient matrices (Hermitian spline)

% ee positions
% index: p{leg}{xyz} = [4 x n] 4 polynomial coeff for n splines 
p = cell(model.NLEGS, 1); 
p_ref = p;
for i = 1:model.NLEGS % n contact points
    for j = 1:3 % xyz indices
        p{i}{j} = opti.variable(4, length(contactState_posns{i})); % x0, x0_dot, x1, x1_dot, 2 splines for each swing phase
        p_ref{i}{j} = ones(4, length(contactState_posns{i}));
    end
end

% ee forces
% index: p{leg}{xyz} = [4 x n] 4 polynomial coeff for n splines
f = cell(model.NLEGS, 1); 
f_ref = f;
for i = 1:model.NLEGS % n contact points
    for j = 1:3 % xyz indices
        f{i}{j} = opti.variable(4, length(contactState_forces{i})); % x0, x0_dot, x1, x1_dot, 3 splines for each contact phase
        f_ref{i}{j} = ones(4, length(contactState_forces{i}));
    end
end

% parameters
r_init      = opti.parameter(3, 1);     % initial CoM position
theta_init  = opti.parameter(3, 1);     % initial CoM orientation
r_des       = opti.parameter(3, 1);     % desired CoM position
theta_des   = opti.parameter(3, 1);     % desired CoM orientation
T           = opti.parameter();         % total time
dt_dyn      = opti.parameter();         % dynamic dt
dt_bs       = opti.parameter();         % base spline duration

% robot/environment parameters
mu          = opti.parameter();          % friction
l_leg_max   = opti.parameter();          % max leg length
f_max       = opti.parameter();          % max force
mass        = opti.parameter();          % mass of robot
Ib          = opti.parameter(3, 1);      % base inertia
Ib_inv      = opti.parameter(3, 1);      % base inertia inverse


%% optimization - generate spline functions

evalSpline_footPositions = cell(model.NLEGS, 3);    % {4 x 3} cell, rows: feet, cols: xyz
evalSpline_footForces = cell(model.NLEGS, 3);       % {4 x 3} cell, rows: feet, cols: xyz

for i = 1:model.NLEGS
    for j = 1:3  % xyz index
        evalSpline_footForces{i, j} = getSplineFunctions(phase_durations{i}, f{i}{j}, ...
                                                         contactDuration_forces{i}, 'Hermite');
        evalSpline_footPositions{i, j} = getSplineFunctions(phase_durations{i}, p{i}{j}, ...
                                                            contactDuration_posns{i}, 'Hermite');
    end
end
evalSpline_comPosn = {getSplineFunctions(phase_durations_base, com_x, [], 'Power')
                      getSplineFunctions(phase_durations_base, com_y, [], 'Power')
                      getSplineFunctions(phase_durations_base, com_z, [], 'Power')};
                  
evalSpline_comOri = {getSplineFunctions(phase_durations_base, roll, [], 'Power')
                     getSplineFunctions(phase_durations_base, pitch, [], 'Power')
                     getSplineFunctions(phase_durations_base, yaw, [], 'Power')};

%% optimization - a-priori constraints from contact sequence
                  
for i = 1:model.NLEGS
    for j = 1:length(contactState_forces{i})                    % force constraints
        if (contactState_forces{i}(j) == 0)                 
            opti.subject_to(f{i}{1}(:, j) == zeros(4, 1));      % in flight, forces are zero
            opti.subject_to(f{i}{2}(:, j) == zeros(4, 1));
            opti.subject_to(f{i}{3}(:, j) == zeros(4, 1));
        end
    end
    for k = 1:length(contactState_posns{i})                     % posn constraints
        if(contactState_forces{i}(k) == 1)
            opti.subject_to(p{i}{1}(2:4, k) == zeros(3, 1));    % foot is constant during contact
            opti.subject_to(p{i}{2}(2:4, k) == zeros(3, 1));  
            opti.subject_to(p{i}{3}(:, k) == zeros(4, 1));      % foot is on ground during contact
        end 
    end
    opti.subject_to(sum(phase_durations{i}) == T);          % phase durations must sum to total time
    opti.subject_to(0 <= phase_durations{i} <= 2.0);               % phase durations must be positive
end

% do I need to link the phase_durations with the spline contact states more
% explicitly? they should be related through the generated functions, but
% not sure how ipopt/casadi will handle it...
                  
%% optimization - constraints

% NOTES:
% terminal constraints might be okay, might not be... double check how they're evaluated

%% initial state constraint
r_0 = [polyval(com_x(:, 1), 0); polyval(com_y(:, 1), 0); polyval(com_z(:, 1), 0)];
theta_0 = [polyval(roll(:, 1), 0); polyval(pitch(:, 1), 0); polyval(yaw(:, 1), 0)];

opti.subject_to(r_0 == r_init);
opti.subject_to(theta_0 == theta_init);
for i = 1:model.NLEGS
    opti.subject_to(evalSpline_footForces{i, 3}.x(dt_dyn*0 - eps) >= 0);
end

%% final state constraints
r_T = [polyval(com_x(:, end), T); polyval(com_y(:, end), T); polyval(com_z(:, end), T)];
theta_T = [polyval(roll(:, end), T); polyval(pitch(:, end), T); polyval(yaw(:, end), T)];

opti.subject_to(r_T >= r_des - 0.02*ones(3,1));
opti.subject_to(r_T <= r_des + 0.02*ones(3,1));
opti.subject_to(theta_T >= theta_des + 0.02*ones(3,1));
opti.subject_to(theta_T <= theta_des + 0.02*ones(3,1));

%% dynamic constraints
q_leg_home = [0 -1.45 2.65];
q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
[~,Ibody_val] = get_mass_matrix(model, q_home, 0);
mass_val = Ibody_val(6,6);
Ibody_inv_val = inv(Ibody_val(1:3,1:3));

for k = 1:N_timesteps-1
    t_k = k*dt_dyn;                                 % current time
    f_k = casadi.MX.zeros(model.NLEGS, 3);          % current foot forces
    p_k = casadi.MX.zeros(model.NLEGS, 3);          % current foot positions
    
    for leg = 1:model.NLEGS
        for xyz = 1:3
            f_k(leg, xyz) = evalSpline_footForces{leg, xyz}.x(t_k);
            p_k(leg, xyz) = evalSpline_footPositions{leg, xyz}.x(t_k);
        end
    end
    
    p_k = reshape(p_k', [model.NLEGS*3, 1]);          % may be incorrect, double check reshaping
    f_k = reshape(f_k', [model.NLEGS*3, 1]);
    
    % non-negative GRF
    opti.subject_to(f_k([3 6 9 12]) >= zeros(4,1));
    opti.subject_to(f_k([3 6 9 12]) <= repmat(f_max,4,1));
    
    % friction Constraints, Eq (7k)
    opti.subject_to(f_k([1 4 7 10]) <= 0.71*mu*f_k([3 6 9 12]));
    opti.subject_to(f_k([1 4 7 10]) >= -0.71*mu*f_k([3 6 9 12]));
    opti.subject_to(f_k([2 5 8 11]) <= 0.71*mu*f_k([3 6 9 12]));
    opti.subject_to(f_k([2 5 8 11]) >= -0.71*mu*f_k([3 6 9 12]));

    r_k         = [evalSpline_comPosn{1}.x(t_k);              % current CoM posn
                   evalSpline_comPosn{2}.x(t_k);
                   evalSpline_comPosn{3}.x(t_k)];
    rDDot_k     = [evalSpline_comPosn{1}.x_ddot(t_k);
                   evalSpline_comPosn{2}.x_ddot(t_k);
                   evalSpline_comPosn{3}.x_ddot(t_k)];
    
    rpy_k       = [evalSpline_comOri{1}.x(t_k);             % current CoM Euler angles
                   evalSpline_comOri{2}.x(t_k);
                   evalSpline_comOri{3}.x(t_k)];
    rpyDot_k    = [evalSpline_comOri{1}.x_dot(t_k);
                   evalSpline_comOri{2}.x_dot(t_k);
                   evalSpline_comOri{3}.x_dot(t_k);];
    rpyDDot_k   = [evalSpline_comOri{1}.x_ddot(t_k);
                   evalSpline_comOri{2}.x_ddot(t_k);
                   evalSpline_comOri{3}.x_ddot(t_k)];
    
    % Appendix B of Winkler RAL paper, angular vel./accel. in world frame
    omega_k = BmatF(rpy_k)*rpyDot_k;
    omegaDot_k = BmatF_dot(rpy_k, rpyDot_k)*rpyDot_k + BmatF(rpy_k)*rpyDDot_k;

    R_world_to_body = rpyToRotMat(rpy_k(1:3))';     % rotation matrices
    R_body_to_world = rpyToRotMat(rpy_k(1:3));
    
    % dynamics
    rddot = (1/mass).*sum(reshape(f_k,3,model.NLEGS),2) + model.gravity';
    omegaDot = diag(Ib_inv)*...
        (R_world_to_body*(cross(p_k(1:3) - r_k(1:3), f_k(1:3))+...
        cross(p_k(4:6) - r_k(1:3), f_k(4:6))+...
        cross(p_k(7:9) - r_k(1:3), f_k(7:9))+...
        cross(p_k(10:12) - r_k(1:3), f_k(10:12)))-...
        cross(R_world_to_body*omega_k(1:3), diag(Ib)*(R_world_to_body*omega_k(1:3))));
    
    opti.subject_to(rDDot_k == rddot);
    % opti.subject_to(omegaDot_k == omegaDot);

%     % kinematic constraints
%     R_yaw = rpyToRotMat([0 0 rpy_k(3)]);
%     for leg = 1:model.NLEGS
%         xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
%         % p_hip = r_k + R_body_to_world*params.hipSrbmLocation(leg,:)';
%         p_hip = r_k + R_yaw*params.hipSrbmLocation(leg,:)';
%         kin_box_dim = 0.05;
%         opti.subject_to((p_k(xyz_idx(1)) - (p_hip(1) + kin_box_dim)) <= 0);
%         opti.subject_to((p_k(xyz_idx(1)) - (p_hip(1) - kin_box_dim)) >= 0);
%         opti.subject_to((p_k(xyz_idx(2)) - (p_hip(2) + kin_box_dim)) <= 0);
%         opti.subject_to((p_k(xyz_idx(2)) - (p_hip(2) - kin_box_dim)) >= 0);
%         % opti.subject_to((p_k(xyz_idx(3)) - (p_hip(3) + 0.375)) >= 0);
%     end
end

%% continuity constraints
t_cont = 0;
for i = 1:N_baseSpline
    t_cont = i*dt_bs;
    for xyz = 1:3
        opti.subject_to(evalSpline_comPosn{xyz}.x(t_cont - eps) == evalSpline_comPosn{xyz}.x(t_cont + eps));
        opti.subject_to(evalSpline_comOri{xyz}.x(t_cont - eps) == evalSpline_comOri{xyz}.x(t_cont + eps));
        opti.subject_to(evalSpline_comPosn{xyz}.x_dot(t_cont - eps) == evalSpline_comPosn{xyz}.x_dot(t_cont + eps));
        opti.subject_to(evalSpline_comOri{xyz}.x_dot(t_cont - eps) == evalSpline_comOri{xyz}.x_dot(t_cont + eps));
        %opti.subject_to(evalSpline_comPosn{xyz}.x_ddot(t_cont - eps) == evalSpline_comPosn{xyz}.x_ddot(t_cont + eps));
        %opti.subject_to(evalSpline_comOri{xyz}.x_ddot(t_cont - eps) == evalSpline_comOri{xyz}.x_ddot(t_cont + eps));
    end
end

for leg = 1:model.NLEGS
    for i = 1:length(contactDuration_forces{leg})
        for xyz = 1:3
            opti.subject_to(evalSpline_footForces{leg, xyz}.x(sum(contactDuration_forces{leg}(1:i)) - eps) == ...
                            evalSpline_footForces{leg, xyz}.x(sum(contactDuration_forces{leg}(1:i)) + eps));
            opti.subject_to(evalSpline_footForces{leg, xyz}.x_dot(sum(contactDuration_forces{leg}(1:i)) - eps) == ...
                            evalSpline_footForces{leg, xyz}.x_dot(sum(contactDuration_forces{leg}(1:i)) + eps));
        end
    end
    for j = 1:length(contactDuration_posns{leg})
        for xyz = 1:3
            opti.subject_to(evalSpline_footPositions{leg, xyz}.x(sum(contactDuration_posns{leg}(1:j)) - eps) == ...
                            evalSpline_footPositions{leg, xyz}.x(sum(contactDuration_posns{leg}(1:j)) + eps));
            opti.subject_to(evalSpline_footPositions{leg, xyz}.x_dot(sum(contactDuration_posns{leg}(1:j)) - eps) == ...
                            evalSpline_footPositions{leg, xyz}.x_dot(sum(contactDuration_posns{leg}(1:j)) + eps));
        end
    end
end

%% optimization - set parameters
theta_init_val = zeros(3,1);
theta_des_val = [0 0 0];
r_init_val = [0 0 0.3]';
r_des_val = [0 0 0.3]';

opti.set_value(r_init, r_init_val);
opti.set_value(theta_init, theta_init_val);
opti.set_value(r_des, r_des_val);
opti.set_value(theta_des, theta_des_val);
opti.set_value(T, T_val);
opti.set_value(dt_dyn, dt_dyn_val);
opti.set_value(dt_bs, dt_baseSpline);

mu_val = 1;
l_leg_max_val = .3;
f_max_val = 200;

opti.set_value(mu, mu_val);
opti.set_value(l_leg_max, l_leg_max_val);
opti.set_value(f_max,f_max_val);
opti.set_value(mass,mass_val);
opti.set_value(Ib,diag(Ibody_val(1:3,1:3)));
opti.set_value(Ib_inv,diag(Ibody_inv_val(1:3,1:3)));

%% set initial guess
c_init_val = repmat(r_init_val(1:3),4,1)+...
    diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -r_init_val(3)],1,4)';
c_init_val = reshape(c_init_val, [3, 4]);

% need to set initial guesses for phase durations to non-zero to not freak
% out IPOPT
for i = 1:model.NLEGS
    N_phases = length(phase_durations{i});
    opti.set_initial(phase_durations{i}, (T_val/N_phases)*ones(1, N_phases));
    for xyz = 1:3
        opti.set_initial(p{i}{xyz}(1, :), c_init_val(xyz, i));
    end
end

% for i = 1:N_baseSpline
%     opti.set_initial(com_x(5,i), r_init_val(1));
%     opti.set_initial(com_y(5,i), r_init_val(2));
%     opti.set_initial(com_z(5,i), r_init_val(3));
% end


%% casadi and IPOPT options
s_opts = struct('max_iter',3000,... %'max_cpu_time',9.0,...
    'tol', 1e-4,... % (1e-6), 1e-4 works well
    'acceptable_tol', 1e-4,... % (1e-4)
    'constr_viol_tol', 1e-3,... % (1e-6), 1e3 works well
    'acceptable_iter', 5,... % (15), % 5 works well
    'nlp_scaling_method','gradient-based',... {'gradient-based','none','equilibration-based'};
    'nlp_scaling_max_gradient',50,... % (100), % 50 works well
    'bound_relax_factor', 1e-6,... % (1e-8), % 1e-6 works well
    'fixed_variable_treatment','relax_bounds',... % {'make_parameter','make_constraint','relax_bounds'}; % relax bounds works well
    'bound_frac',5e-1,... % (1e-2), 5e-1 works well
    'bound_push',5e-1,... % (1e-2), 5e-1 works well
    'mu_strategy','adaptive',... % {'monotone','adaptive'}; % adaptive works very well
    'mu_oracle','probing',... % {'quality-function','probing','loqo'}; % probing works very well
    'fixed_mu_oracle','probing',... % {'average_compl','quality-function','probing','loqo'}; % probing decent
    'adaptive_mu_globalization','obj-constr-filter',... % {'obj-constr-filter','kkt-error','never-monotone-mode'};
    'mu_init',1e-1,... % [1e-1 1e-2 1]
    'alpha_for_y','bound-mult',... % {'primal','bound-mult','min','max','full','min-dual-infeas','safer-min-dual-infeas','primal-and-full'}; % primal or bound-mult seems best
    'alpha_for_y_tol',1e1,... % (1e1)
    'recalc_y','no',... % {'no','yes'};
    'max_soc',4,... % (4)
    'accept_every_trial_step','no',... % {'no','yes'}
    'linear_solver','mumps',... % {'ma27','mumps','ma57','ma77','ma86'} % ma57 seems to work well
    'linear_system_scaling','slack-based',... {'mc19','none','slack-based'}; % Slack-based
    'linear_scaling_on_demand','yes',... % {'yes','no'};
    'max_refinement_steps',10,... % (10)
    'min_refinement_steps',1,... % (1)
    'warm_start_init_point', 'no'); % (no)

s_opts.file_print_level = 2;
s_opts.print_level = 4;
s_opts.print_frequency_iter = 100;
s_opts.print_timing_statistics ='no';

opti.solver('ipopt', struct(), s_opts);


%% solve
disp_box('Solving with Opti Stack');
sol = opti.solve_limited();
% opti.debug.show_infeasibilities(1e-3)

%% unpacking solution

p_star = cell(model.NLEGS, 1); 
f_star = cell(model.NLEGS, 1); 
phaseDurations_star = cell(model.NLEGS, 1);
contactDurations_force_star = cell(model.NLEGS, 1);
contactDurations_posn_star = cell(model.NLEGS, 1);
for i = 1:model.NLEGS
    for xyz = 1:3
        p_star{i}{xyz} = sol.value(p{i}{xyz});
        f_star{i}{xyz} = sol.value(f{i}{xyz});
    end
    phaseDurations_star{i} = sol.value(phase_durations{i});
    contactDurations_force_star{i} = sol.value(contactDuration_forces{i});
    contactDurations_posn_star{i} = sol.value(contactDuration_posns{i});
end

q_star = cell(6, 1);
q_star{1} = sol.value(com_x); q_star{2} = sol.value(com_y); q_star{3} = sol.value(com_z);
q_star{4} = sol.value(roll); q_star{5} = sol.value(pitch); q_star{6} = sol.value(yaw);

%% evaluation 

t_star = 0:0.01:T_val;
q_star_eval = zeros(6, length(t_star));
f_star_eval = cell(model.NLEGS, 1); 

for i = 1:length(t_star)
    for j = 1:6
        q_star_eval(j,i) = evaluateSpline(q_star{j}, phase_durations_base, t_star(i), 'Power');
    end
    for legs = 1:model.NLEGS
        for xyz = 1:3
            f_star_eval{legs}(xyz, i) = evaluateSpline(f_star{legs}{xyz}, contactDurations_force_star{legs}, t_star(i), 'Hermite');
            p_star_eval{legs}(xyz, i) = evaluateSpline(p_star{legs}{xyz}, contactDurations_posn_star{legs}, t_star(i), 'Hermite');
        end
    end
end

%% visualization and plots

figure;
hold on;
plot(t_star, q_star_eval(1,:), 'r--')
plot(t_star, q_star_eval(2,:), 'g--')
plot(t_star, q_star_eval(3,:), 'b--')
legend('x', 'y', 'z')
xlabel('Time (s)')
ylabel('Posn (m)')

figure;
hold on;
plot(t_star, f_star_eval{1}(3, :), 'r--')
plot(t_star, f_star_eval{2}(3, :), 'g--')
plot(t_star, f_star_eval{3}(3, :), 'b--')
plot(t_star, f_star_eval{4}(3, :), 'k--')
legend('FR', 'FL', 'BR', 'BL')
xlabel('Time (s)')
ylabel('F_z (N)')

figure;
hold on;
plot(t_star, f_star_eval{1}(1, :), 'r--')
plot(t_star, f_star_eval{2}(1, :), 'g--')
plot(t_star, f_star_eval{3}(1, :), 'b--')
plot(t_star, f_star_eval{4}(1, :), 'k--')
legend('FR', 'FL', 'BR', 'BL')
xlabel('Time (s)')
ylabel('F_x (N)')

figure;
hold on;
plot(t_star, f_star_eval{1}(2, :), 'r--')
plot(t_star, f_star_eval{2}(2, :), 'g--')
plot(t_star, f_star_eval{3}(2, :), 'b--')
plot(t_star, f_star_eval{4}(2, :), 'k--')
legend('FR', 'FL', 'BR', 'BL')
xlabel('Time (s)')
ylabel('F_y (N)')

figure;
hold on;
plot(t_star, p_star_eval{1}(3, :), 'r--')
plot(t_star, p_star_eval{2}(3, :), 'g--')
plot(t_star, p_star_eval{3}(3, :), 'b--')
plot(t_star, p_star_eval{4}(3, :), 'k--')
legend('FR', 'FL', 'BR', 'BL')
xlabel('Time (s)')
ylabel('Foot height (m)')

figure;
hold on;
plot(t_star, p_star_eval{1}(1, :), 'r--')
plot(t_star, p_star_eval{2}(1, :), 'g--')
plot(t_star, p_star_eval{3}(1, :), 'b--')
plot(t_star, p_star_eval{4}(1, :), 'k--')
legend('FR', 'FL', 'BR', 'BL')
xlabel('Time (s)')
ylabel('Foot xPosn (m)')


fb_motion = [q_star_eval; repmat(repmat(q_leg_home', 4, 1),1,length(t_star))];
if show_animation
    showmotion(model,t_star,fb_motion)
end


