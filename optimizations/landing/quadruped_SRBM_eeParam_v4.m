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
addpath(genpath('../../utilities_general'));
addpath(genpath('utilities_eeParam'));
import casadi.*

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModelMC3D(params, model, 0);

%% problem parameters and setup

dt_dyn = 0.1;                                               % timestep for dynamics
T_val = 1.0;                                                % total time
N_timesteps = round(T_val/dt_dyn + 1);                      % number of steps for dynamics

dt_base = 0.2;                                              % duration of polynomials for base motion
N_base = round(T_val/dt_base);                              % number of polynomials for base
time_base = zeros(1, N_base + 1);                           % time breakpoints of base piecewise polynomial
order_base = 4;                                             % polynomial order of base splines
for i = 1: N_base
    time_base(i + 1) = dt_base*i;
end

% contact sequence specification
% [# contact phases, # flight phases, begin_in_contact (1 or 0)]
contactPhases = [1 1 0;         % FR foot  
                 1 1 1;         % FL foot
                 1 1 1;         % BR foot
                 1 1 0];        % BL foot
contactSequence = generateContactSequence(contactPhases);   % contact sequence for legs
contactState_forces = cell(model.NLEGS, 1);                 % contact state for force splines
contactState_posns = cell(model.NLEGS, 1);                  % contact state for posn splines
contactDuration_forces = contactState_forces;               % contact durations for force splines
contactDuration_posns = contactState_posns;                 % contact durations for posn splines

N_splines_stance = 3;                                       % number of forces splines per stance phase 
N_splines_swing = 2;                                        % number of posn splines per swing phase

%% optimization - variables

opti = casadi.Opti();

% decision variables:

% CoM states, 3rd order polynomial parametrization, fixed duration
% struct: each 'node' contains coef_'state' [xyz, coef] at node i
nodes_base = cell(1, N_base);
for i = 1:N_base
    nodes_base{i}.coef_lin = opti.variable(3, order_base);
    nodes_base{i}.coef_ang = opti.variable(3, order_base);
    nodes_base{i}.coef_accel = getDerivCoef(getDerivCoef(nodes_base{i}.coef_lin));
    nodes_base{i}.coef_alpha = getDerivCoef(getDerivCoef(nodes_base{i}.coef_ang));
    nodes_base{i}.dt = dt_base;
    nodes_base{i}.idx = i;
end

% contact durations
phaseDurations = cell(model.NLEGS, 1);
for i = 1:model.NLEGS
    phaseDurations{i} = opti.variable(1, numel(contactSequence{i}));
end

% contact states and durations of splines
% can reduce # of variables by only defining a single constant for stance/swing phases for posn/force, but keeping for now
time_force = cell(model.NLEGS, 1);
time_posn = cell(model.NLEGS, 1);
for i = 1:model.NLEGS
    temp_forces = [];
    temp_posns = [];
    temp_f_dur = [];
    temp_p_dur = [];
    for j = 1:length(contactSequence{i})
        if (contactSequence{i}(j) == 0)
            % not in contact
            temp_forces = [temp_forces 0];
            temp_posns = [temp_posns zeros(1, N_splines_stance)];
            temp_f_dur = [temp_f_dur phaseDurations{i}(j)];
            temp_p_dur = [temp_p_dur ones(1,N_splines_swing)*(phaseDurations{i}(j)/N_splines_swing)];
        else
            % in contact
            temp_forces = [temp_forces ones(1, 3)];
            temp_posns = [temp_posns 1];
            temp_f_dur = [temp_f_dur ones(1,N_splines_stance)*(phaseDurations{i}(j)/N_splines_stance)];
            temp_p_dur = [temp_p_dur phaseDurations{i}(j)];
        end
    end
    contactState_forces{i} = temp_forces;
    contactState_posns{i} = temp_posns;
    contactDuration_forces{i} = temp_f_dur;
    contactDuration_posns{i} = temp_p_dur;
    time_force{i} = casadi.MX.zeros(1, numel(contactState_forces{i}) + 1);
    time_posn{i} = casadi.MX.zeros(1, numel(contactState_posns{i}) + 1);
end

nodes_force = cell(model.NLEGS, 1);
nodes_posn = cell(model.NLEGS, 1);
for leg = 1:model.NLEGS
    nodes_force{leg}.x = opti.variable(numel(contactState_forces{leg}), 4);
    nodes_force{leg}.y = opti.variable(numel(contactState_forces{leg}), 4);
    nodes_force{leg}.z = opti.variable(numel(contactState_forces{leg}), 4);
    nodes_force{leg}.duration = contactDuration_forces{leg};
    nodes_force{leg}.cs = contactState_forces{leg};
    
    nodes_posn{leg}.x = opti.variable(numel(contactState_posns{leg}), 4);
    nodes_posn{leg}.y = opti.variable(numel(contactState_posns{leg}), 4);
    nodes_posn{leg}.z = opti.variable(numel(contactState_posns{leg}), 4);
    nodes_posn{leg}.duration = contactDuration_forces{leg};
    nodes_posn{leg}.cs = contactState_posns{leg};
    
    for i = 1:numel(contactState_forces{leg})
        time_force{leg}(1, i + 1) = sum(contactDuration_forces{leg}(1:i));
    end
    for j = 1:numel(contactState_posns{leg})
        time_posn{leg}(1, j + 1) = sum(contactDuration_forces{leg}(1:j));
    end
    
end

% initial + terminal parameters
r_init          = opti.parameter(3, 1);     % initial CoM position
r_des           = opti.parameter(3, 1);     % desired CoM position
rDot_init       = opti.parameter(3, 1);     % initial CoM l. velocity
rDot_des        = opti.parameter(3, 1);     % desired CoM l. velocity
theta_init      = opti.parameter(3, 1);     % initial CoM orientation
theta_des       = opti.parameter(3, 1);     % desired CoM orientation
thetaDot_init   = opti.parameter(3, 1);     % initial CoM a. velocity
thetaDot_des    = opti.parameter(3, 1);     % desired CoM a. velocity

% time parameters
T               = opti.parameter();         % total time
dt_test         = opti.parameter();
opti.set_value(dt_test, 0.1);
% check if needed
% dt_dyn          = opti.parameter();         % dynamic dt  
% dt_bs           = opti.parameter();         % base spline duration

% robot/environment parameters
mu          = opti.parameter();          % friction
l_leg_max   = opti.parameter();          % max leg length
f_max       = opti.parameter();          % max force
mass        = opti.parameter();          % mass of robot
Ib          = opti.parameter(3, 1);      % base inertia
Ib_inv      = opti.parameter(3, 1);      % base inertia inverse

%% optimization - a-priori constraints from contact sequence

for leg = 1:model.NLEGS
    for i = 1:numel(nodes_force{leg}.cs)
        f_0 = [nodes_force{leg}.x(i, 1); nodes_force{leg}.y(i, 1); nodes_force{leg}.z(i, 1)]; % x0 of xyz force spline
        f_1 = [nodes_force{leg}.x(i, 3); nodes_force{leg}.y(i, 3); nodes_force{leg}.z(i, 3)]; % x1 of xyz force spline
        if (nodes_force{leg}.cs(i) == 0)
            % alternatively, could constrain f_i == 0? unsure which is more efficient
            opti.subject_to(nodes_force{leg}.x(i, :) == zeros(1, 4));
            opti.subject_to(nodes_force{leg}.y(i, :) == zeros(1, 4));
            opti.subject_to(nodes_force{leg}.z(i, :) == zeros(1, 4));
        else
            opti.subject_to(f_max >= f_0(3) >= 0);                                 % ground normals >= 0 in contact
            opti.subject_to(-0.71*f_0(3) <= f_0(1) <= 0.71*mu*f_0(3));    % friction pyramid x
            opti.subject_to(-0.71*f_0(3) <= f_0(2) <= 0.71*mu*f_0(3));    % friction pyramid y
            opti.subject_to(f_max >= f_1(3) >= zeros(3, 1));   % needed? with continuity, seems unnecessary...
        end
    end
    for j = 1:numel(nodes_posn{leg}.cs)
        if (nodes_posn{leg}.cs(j) == 1)
            opti.subject_to(nodes_posn{leg}.x(j,1) == nodes_posn{leg}.x(j,3));
            opti.subject_to(nodes_posn{leg}.y(j,1) == nodes_posn{leg}.y(j,3));
            opti.subject_to(nodes_posn{leg}.z(j,:) == zeros(1, 4));
        else
            opti.subject_to(nodes_posn{leg}.z(j,:) >= zeros(1, 4));
        end        
    end
    opti.subject_to(sum(phaseDurations{leg}) == T);               % phase durations of legs sum to total time
    for k = 1:numel(contactSequence{leg})
        opti.subject_to(0.25 <= phaseDurations{leg}(k) <= T);        % phase durations must be positive and bounded
    end
end

%% optimization - constraints

%% initial state constraint
r_0     = [polyval(nodes_base{1}.coef_lin(1, :)', 0);
           polyval(nodes_base{1}.coef_lin(2, :)', 0);
           polyval(nodes_base{1}.coef_lin(3, :)', 0)];
r_0_dot = [polyval(getDerivCoef(nodes_base{1}.coef_lin(1, :))', 0);
           polyval(getDerivCoef(nodes_base{1}.coef_lin(2, :))', 0);
           polyval(getDerivCoef(nodes_base{1}.coef_lin(3, :))', 0)];
theta_0 = [polyval(nodes_base{1}.coef_ang(1, :)', 0);
           polyval(nodes_base{1}.coef_ang(2, :)', 0);
           polyval(nodes_base{1}.coef_ang(3, :)', 0)];
theta_0_dot = [polyval(getDerivCoef(nodes_base{1}.coef_ang(1, :))', 0);
           polyval(getDerivCoef(nodes_base{1}.coef_ang(2, :))', 0);
           polyval(getDerivCoef(nodes_base{1}.coef_ang(3, :))', 0)];
r_0_accel     = [polyval(nodes_base{1}.coef_accel(1, :)', 0);
           polyval(nodes_base{1}.coef_accel(2, :)', 0);
           polyval(nodes_base{1}.coef_accel(3, :)', 0)];

opti.subject_to(r_0 == r_init);
opti.subject_to(r_0_dot == rDot_init);
opti.subject_to(theta_0 == theta_init);
opti.subject_to(theta_0_dot == thetaDot_init);

opti.subject_to(r_0_accel == zeros(3,1));

%% final state constraints
r_f     = [polyval(nodes_base{end}.coef_lin(1, :)', dt_base);
           polyval(nodes_base{end}.coef_lin(2, :)', dt_base);
           polyval(nodes_base{end}.coef_lin(3, :)', dt_base)];
r_f_dot     = [polyval(getDerivCoef(nodes_base{end}.coef_lin(1, :))', dt_base);
           polyval(getDerivCoef(nodes_base{end}.coef_lin(2, :))', dt_base);
           polyval(getDerivCoef(nodes_base{end}.coef_lin(3, :))', dt_base)];       
theta_f = [polyval(nodes_base{end}.coef_ang(1, :)', dt_base);
           polyval(nodes_base{end}.coef_ang(2, :)', dt_base);
           polyval(nodes_base{end}.coef_ang(3, :)', dt_base)];

opti.subject_to(r_f >= r_des - 0.02*ones(3,1));
opti.subject_to(r_f <= r_des + 0.02*ones(3,1));
opti.subject_to(theta_f >= theta_des + 0.02*ones(3,1));
opti.subject_to(theta_f <= theta_des + 0.02*ones(3,1));

opti.subject_to(r_f_dot == zeros(3, 1));

%% continuity constraints
for i = 1:numel(nodes_base)-1
    for xyz = 1:3
        prevSpline_end_lin = polyval(nodes_base{i}.coef_lin(xyz, :)', dt_base);
        prevSpline_end_ang = polyval(nodes_base{i}.coef_ang(xyz, :)', dt_base);
        prevSpline_end_v = polyval(getDerivCoef(nodes_base{i}.coef_lin(xyz, :))', dt_base);
        prevSpline_end_w = polyval(getDerivCoef(nodes_base{i}.coef_lin(xyz, :))', dt_base);
        prevSpline_end_accel = polyval(nodes_base{i}.coef_accel(xyz, :)', dt_base);
        prevSpline_end_alpha = polyval(nodes_base{i}.coef_alpha(xyz, :)', dt_base);
        nextSpline_start_lin = polyval(nodes_base{i+1}.coef_lin(xyz, :)', 0);
        nextSpline_start_ang = polyval(nodes_base{i+1}.coef_ang(xyz, :)', 0);
        nextSpline_start_v = polyval(getDerivCoef(nodes_base{i+1}.coef_lin(xyz, :))', 0);
        nextSpline_start_w = polyval(getDerivCoef(nodes_base{i+1}.coef_ang(xyz, :))', 0);
        nextSpline_start_accel = polyval(nodes_base{i+1}.coef_accel(xyz, :)', 0);
        nextSpline_start_alpha = polyval(nodes_base{i+1}.coef_alpha(xyz, :)', 0);
        
        opti.subject_to(prevSpline_end_lin == nextSpline_start_lin);
        opti.subject_to(prevSpline_end_ang == nextSpline_start_ang);
        opti.subject_to(prevSpline_end_v == nextSpline_start_v);
        opti.subject_to(prevSpline_end_w == nextSpline_start_w);
        opti.subject_to(prevSpline_end_accel == nextSpline_start_accel);
        opti.subject_to(prevSpline_end_alpha == nextSpline_start_alpha);
    end
end
% dont think we need to constraint velocities? spline 'smoothness' likely
% comes from dynamics? but cont. accel doesnt guarantee cont. velocities...

for leg = 1:model.NLEGS
    for i = 1:numel(nodes_force{leg}.cs)-1
        x1_prev = [nodes_force{leg}.x(i, 3); nodes_force{leg}.y(i, 3); nodes_force{leg}.z(i, 3)];                 % x1 of prev force spline
        x0_next = [nodes_force{leg}.x(i+1, 1); nodes_force{leg}.y(i+1, 1); nodes_force{leg}.z(i+1, 1)];           % x0 of next force spline
        x1_dot_prev = [nodes_force{leg}.x(i, 4); nodes_force{leg}.y(i, 4); nodes_force{leg}.z(i, 4)];             % x1_dot of prev force spline
        x0_dot_next = [nodes_force{leg}.x(i, 2); nodes_force{leg}.y(i, 2); nodes_force{leg}.z(i, 2)];             % x0_dot of next force spline
        
        opti.subject_to(x1_prev == x0_next);
        opti.subject_to(x1_dot_prev == x0_dot_next);
    end
    for j = 1:numel(nodes_posn{leg}.cs)-1
        x1_prev = [nodes_posn{leg}.x(j, 3); nodes_posn{leg}.y(j, 3); nodes_posn{leg}.z(j, 3)];                    % x1 of prev posn spline
        x0_next = [nodes_posn{leg}.x(j+1, 1); nodes_posn{leg}.y(j+1, 1); nodes_posn{leg}.z(j+1, 1)];              % x0 of next posn spline
        x1_dot_prev = [nodes_posn{leg}.x(j, 4); nodes_posn{leg}.y(j, 4); nodes_posn{leg}.z(j, 4)];                % x1_dot of prev posn spline
        x0_dot_next = [nodes_force{leg}.x(i, 2); nodes_force{leg}.y(i, 2); nodes_force{leg}.z(i, 2)];             % x0_dot of next posn spline
        
        opti.subject_to(x1_dot_prev == x0_dot_next);
        opti.subject_to(x1_prev == x0_next);
    end
end


%% dynamic constraints
q_leg_home = [0 -1.45 2.65];
q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
[~,Ibody_val] = get_mass_matrix(model, q_home, 0);
mass_val = Ibody_val(6,6);
Ibody_inv_val = inv(Ibody_val(1:3,1:3));

for k = 1:N_timesteps-1
    t_global = k*dt_dyn;                                 % current time
    t_global_test = k*dt_test;
    base_idx = full(getSplineIndex(time_base, t_global));
    t_local = t_global - time_base(base_idx);

    f_k = casadi.MX.zeros(3*model.NLEGS, 1);
    p_k = casadi.MX.zeros(3*model.NLEGS, 1);

    for leg = 1:model.NLEGS
        f_idx = full(getSplineIndex(time_force{leg}, t_global_test));
        p_idx = full(getSplineIndex(time_posn{leg}, t_global_test));
        
        t_local_f = t_global - time_force{leg}(f_idx);
        f_k(3*(leg-1) + 1) = polyval(convertHermiteCoef(nodes_force{leg}.x(f_idx, :), nodes_force{leg}.duration(f_idx))', t_local_f);
        f_k(3*(leg-1) + 2) = polyval(convertHermiteCoef(nodes_force{leg}.y(f_idx, :), nodes_force{leg}.duration(f_idx))', t_local_f);
        f_k(3*(leg-1) + 3) = polyval(convertHermiteCoef(nodes_force{leg}.z(f_idx, :), nodes_force{leg}.duration(f_idx))', t_local_f);
        
        t_local_p = t_global - time_posn{leg}(p_idx);
        p_k(3*(leg-1) + 1) = polyval(convertHermiteCoef(nodes_posn{leg}.x(p_idx, :), nodes_posn{leg}.duration(p_idx))', t_local_p);
        p_k(3*(leg-1) + 2) = polyval(convertHermiteCoef(nodes_posn{leg}.y(p_idx, :), nodes_posn{leg}.duration(p_idx))', t_local_p);
        p_k(3*(leg-1) + 3) = polyval(convertHermiteCoef(nodes_posn{leg}.z(p_idx, :), nodes_posn{leg}.duration(p_idx))', t_local_p);
    end

    
    r_k     = [polyval(nodes_base{base_idx}.coef_lin(1, :)', t_local);
               polyval(nodes_base{base_idx}.coef_lin(2, :)', t_local);
               polyval(nodes_base{base_idx}.coef_lin(3, :)', t_local)];
    rDDot_k = [polyval(nodes_base{base_idx}.coef_accel(1, :)', t_local);
               polyval(nodes_base{base_idx}.coef_accel(2, :)', t_local);
               polyval(nodes_base{base_idx}.coef_accel(3, :)', t_local)];
           
    rpy_k       = [polyval(nodes_base{base_idx}.coef_ang(1, :)', t_local);
                   polyval(nodes_base{base_idx}.coef_ang(2, :)', t_local);
                   polyval(nodes_base{base_idx}.coef_ang(3, :)', t_local)];
    rpyDot_k    = [polyval(getDerivCoef(nodes_base{base_idx}.coef_ang(1, :))', t_local);
                   polyval(getDerivCoef(nodes_base{base_idx}.coef_ang(2, :))', t_local);
                   polyval(getDerivCoef(nodes_base{base_idx}.coef_ang(3, :))', t_local)];
    rpyDDot_k   = [polyval(nodes_base{base_idx}.coef_alpha(1, :)', t_local);
                   polyval(nodes_base{base_idx}.coef_alpha(2, :)', t_local);
                   polyval(nodes_base{base_idx}.coef_alpha(3, :)', t_local)];
    
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
    opti.subject_to(omegaDot_k == omegaDot);
    
    for leg = 1:model.NLEGS
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        
        % kinematic Limits - applied only at touchdown
        % do these account for leg lenghts? or is that left to the timing
        % of the contact schedule?
        
        r_hip = r_k(1:3) + R_body_to_world*params.hipSrbmLocation(leg,:)';
        p_rel = R_world_to_body*(p_k(xyz_idx) - r_hip);
        kin_box_x = 0.05;
        kin_box_y = 0.05;
        kin_box_z = 0.30;
        
        opti.subject_to(-kin_box_x <= p_rel(1) <= kin_box_x);
        opti.subject_to(-kin_box_y <= p_rel(2) <= kin_box_y);
        opti.subject_to(-kin_box_z <= p_rel(3) + 0.05 <= 0)
        
    end

end

%% optimization - set parameters
theta_init_val = zeros(3,1);
thetaDot_init_val = zeros(3,1);
theta_des_val = [0 0 0];
r_init_val = [0.5 0 0.3]';
rDot_init_val = [0 0 0]';
r_des_val = [0 0 0.3]';

opti.set_value(r_init, r_init_val);
opti.set_value(rDot_init, rDot_init_val);
opti.set_value(r_des, r_des_val);
opti.set_value(theta_init, theta_init_val);
opti.set_value(thetaDot_init, thetaDot_init_val);
opti.set_value(theta_des, theta_des_val);
opti.set_value(T, T_val);
%opti.set_value(dt_dyn, dt_dyn_val);
%opti.set_value(dt_bs, dt_baseSpline);

mu_val = 1;
l_leg_max_val = .35;
f_max_val = 200;

opti.set_value(mu, mu_val);
opti.set_value(l_leg_max, l_leg_max_val);
opti.set_value(f_max,f_max_val);
opti.set_value(mass,mass_val);
opti.set_value(Ib,diag(Ibody_val(1:3,1:3)));
opti.set_value(Ib_inv,diag(Ibody_inv_val(1:3,1:3)));

%% set initial guess
% need to set initial guesses for phase durations to non-zero to not freak
% out IPOPT
for leg = 1:model.NLEGS
    N_phases = length(phaseDurations{leg});
    opti.set_initial(phaseDurations{leg}, (T_val/N_phases)*ones(1, N_phases));
    
    for i = 1:numel(nodes_force{leg}.cs)
        if nodes_force{leg}.cs(i) == 1
            opti.set_initial(nodes_force{leg}.z(i, 1), 20);
        end
    end
    for j = 1:numel(nodes_posn{leg}.cs)
            opti.set_initial(nodes_posn{leg}.z(j, 1), 0.);

    end
end

cost = casadi.MX(0);
for leg=1:model.NLEGS
    for i = 1:numel(nodes_force{leg}.cs)
        coeff = dot(nodes_force{leg}.x(i, :), nodes_force{leg}.x(i, :));
        cost = cost + coeff^2;
    end
end
opti.minimize(cost);

%% casadi and IPOPT options
p_opts = struct('expand',true); % this speeds up ~x10
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

%% unpack solution

baseCoef_star = cell(6, 1);
for i = 1:numel(baseCoef_star)
    baseCoef_star{i} = zeros(numel(nodes_base), order_base);
end

for i = 1:numel(nodes_base)
    baseCoef_star{1}(i, :) = sol.value(nodes_base{i}.coef_lin(1, :));
    baseCoef_star{2}(i, :) = sol.value(nodes_base{i}.coef_lin(2, :));
    baseCoef_star{3}(i, :) = sol.value(nodes_base{i}.coef_lin(3, :));
    baseCoef_star{4}(i, :) = sol.value(nodes_base{i}.coef_ang(1, :));
    baseCoef_star{5}(i, :) = sol.value(nodes_base{i}.coef_ang(2, :));
    baseCoef_star{6}(i, :) = sol.value(nodes_base{i}.coef_ang(3, :));
end


f_star = cell(3, model.NLEGS);
p_star = cell(3, model.NLEGS);
dur_f = cell(model.NLEGS, 1);
dur_p = cell(model.NLEGS, 1);

for leg = 1:model.NLEGS
    for i = 1:numel(nodes_force{leg}.cs)
        f_star{1, leg}(i, :) = convertHermiteCoef(sol.value(nodes_force{leg}.x(i, :)), sol.value(nodes_force{leg}.duration(i)));
        f_star{2, leg}(i, :) = convertHermiteCoef(sol.value(nodes_force{leg}.y(i, :)), sol.value(nodes_force{leg}.duration(i)));
        f_star{3, leg}(i, :) = convertHermiteCoef(sol.value(nodes_force{leg}.z(i, :)), sol.value(nodes_force{leg}.duration(i)));
    end
    for j = 1:numel(nodes_posn{leg}.cs)
        p_star{1, leg}(j, :) = convertHermiteCoef(sol.value(nodes_posn{leg}.x(j, :)), sol.value(nodes_posn{leg}.duration(j)));
        p_star{2, leg}(j, :) = convertHermiteCoef(sol.value(nodes_posn{leg}.y(j, :)), sol.value(nodes_posn{leg}.duration(j)));
        p_star{3, leg}(j, :) = convertHermiteCoef(sol.value(nodes_posn{leg}.z(j, :)), sol.value(nodes_posn{leg}.duration(j)));
    end
    dur_f{leg} = sol.value(nodes_force{leg}.duration);
    dur_p{leg} = sol.value(nodes_posn{leg}.duration);
end


%% plots
pp = mkpp(time_base, baseCoef_star{3});
time_eval = 0:0.01:T_val;
figure;
plot(time_eval, ppval(pp, time_eval))

figure;
hold on;
pp_f = mkpp(sol.value(time_force{2}), f_star{1, 2});
plot(time_eval, ppval(pp_f, time_eval))
pp_f = mkpp(sol.value(time_force{2}), f_star{2, 2});
plot(time_eval, ppval(pp_f, time_eval))
pp_f = mkpp(sol.value(time_force{2}), f_star{3, 2});
plot(time_eval, ppval(pp_f, time_eval))
xlabel('Time (s)')
ylabel('Force FL (N)')
legend('F_x', 'F_y', 'F_z')
hold off

figure;
pp_p = mkpp(sol.value(time_posn{2}), p_star{3, 2});
plot(time_eval, ppval(pp_p, time_eval));

%% visualization

% fb_motion = [q_star_eval; repmat(repmat(q_leg_home', 4, 1),1,length(t_star))];
% if show_animation
%     showmotion(model,t_star,fb_motion)
% end


