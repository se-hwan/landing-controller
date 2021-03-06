%% metadata
% Description:  Trajectory optimization for quadrupedal landing with single rigid body model
%               Uses Michael Posa's contact complementarity constraints with no-slip contacts
% Author:       Se Hwan Jeon

%% cleanup
rmpath(genpath('.')); % clear all previously added paths
clear; clc; close all;

%% flags
show_animation = true;
run_IK = true;
make_plots = false;

%% add library paths
% may need to specify os directory
% addpath(genpath('../../utilities/casadi/casadi_windows'));
% addpath(genpath('../../utilities/casadi/casadi_linux'));
addpath(genpath('../../../utilities_general'));
addpath(genpath('../codegen_casadi'));
import casadi.*

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModelMC3D(params, model, 0);

%% contact schedule parameters
N = 21; 
T = 0.6;
dt_val = repmat(T/(N-1),1,N-1);

%% optimization
% Optimization variables
% q       (6) [x,y,z,roll,pitch,yaw] (floating base state)
% qdot    (6) [omega (body frame), vBody (world frame)]
% c       (NLEGS*3) (absolute world frame)
% f_grf   (NLEGS*3) (absolute world frame)

opti = casadi.Opti();
X = opti.variable(12, N);               % floating base
q       = X(1:6,:);
qdot    = X(7:12,:);
U = opti.variable(6*model.NLEGS, N-1);  % foot posns + GRFs
c     = U(1:12,:);
f_grf = U(13:24,:);

% Optimization Parameters
Xref = opti.parameter(12, N);       % floating base reference
Uref = opti.parameter(24, N-1);     % foot posn + GRF reference

dt = opti.parameter(1, N-1);                % timesteps

q_min = opti.parameter(6,1);q_max = opti.parameter(6,1);
qd_min = opti.parameter(6,1);qd_max = opti.parameter(6,1);

q_init = opti.parameter(6,1);
qd_init = opti.parameter(6,1);
c_init = opti.parameter(3*model.NLEGS,1);

q_term_min = opti.parameter(6,1);q_term_max = opti.parameter(6,1);
qd_term_min = opti.parameter(6,1);qd_term_max = opti.parameter(6,1);
QX = opti.parameter(12,1);QN = opti.parameter(12,1);    % weighting matrices
Qc = opti.parameter(3,1);Qf = opti.parameter(3,1);

mu = opti.parameter();          % robot/environment parameters
l_leg_max = opti.parameter();
f_max = opti.parameter();
mass = opti.parameter();
Ib = opti.parameter(3,1);
Ib_inv = opti.parameter(3,1);

p_hip = [0.19;-0.1;-0.2;...
         0.19;0.1;-0.2;...
         -0.19;-0.1;-0.2;...
         -0.19;0.1;-0.2]; % TODO: make this an opt parameter?
     
%% cost function
cost = casadi.MX(0);             % initialize cost
for k = 1:(N-1)                  % running cost
%     X_err = X(:,k) - Xref(:,k);                                         % floating base error
%     pf_err = repmat(X(1:3,k),model.N_GND_CONTACTS,1) + p_hip - c(:,k);  % foot position error
%     U_err = U(13:24,k) - Uref(13:24,k);                                 % GRF error
%     cost = cost + (X_err'*diag(QX)*X_err+...                            % sum of quadratic error costs
%         pf_err'*diag(repmat(Qc,4,1))*pf_err+...
%         U_err'*diag(repmat(Qf,4,1))*U_err)*dt(k);
%     
%     if (k < 10)
%         v_world = X(10:12, k);
%         phi_v_world = atan(v_world(1)/v_world(3));
%         R_body_to_world = rpyToRotMat(X(4:6, k));
%         for leg = 1:2
%             xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
%             r_hip = X(1:3)' + R_body_to_world*params.hipSrbmLocation(leg,:)';
%             p_foot_rel = U(xyz_idx, k) - r_hip;
%             phi_p_foot = atan(p_foot_rel(1)/p_foot_rel(3));
%             % cost = cost + .1*(phi_v_world - phi_p_foot)^2;
%         end
%     end
%     U_err = U(13:24,k) - Uref(13:24,k);                                 % GRF error
%     cost = cost +  U_err'*diag(repmat(Qf,4,1))*U_err*dt(k);
end
X_err = X(:,end)-Xref(:,end);    % terminal cost
cost = cost + X_err'*diag(QN)*X_err;

opti.minimize(cost);             % set objective

%% initial state constraint
opti.subject_to(q(1:6,1) == q_init);        % initial pos + ori
opti.subject_to(qdot(1:6,1) == qd_init);    % initial ang. vel. + lin. vel.
% opti.subject_to(c(:,1) == c_init);          % initial foot positions

%% terminal state constraints
opti.subject_to(q(:,N) >= q_term_min);      % bounds terminal state to be within specified min/max values
opti.subject_to(q(:,N) <= q_term_max);
opti.subject_to(qdot(:,N) >= qd_term_min);
opti.subject_to(qdot(:,N) <= qd_term_max);

%% general constraints
q_leg_home = [0 -1.45 2.65];
q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
[~,Ibody_val] = get_mass_matrix(model, q_home, 0);
mass_val = Ibody_val(6,6);
Ibody_inv_val = inv(Ibody_val(1:3,1:3));


for k = 1:N-1               % the 'k' suffix indicates the value of the variable at the current timestep
    
    qk = q(:,k);
    qdk = qdot(:,k);
    rpyk = q(4:6,k);
    ck = c(:,k);
    fk = f_grf(:,k);
    
    R_world_to_body = rpyToRotMat(rpyk(1:3))';
    R_body_to_world = rpyToRotMat(rpyk(1:3));
    
    % dynamics
    rddot = (1/mass).*sum(reshape(fk,3,model.NLEGS),2)+model.gravity';
    omegaDot = diag(Ib_inv)*...
        (R_world_to_body*(cross(ck(1:3)-qk(1:3),fk(1:3))+...
        cross(ck(4:6)-qk(1:3),fk(4:6))+...
        cross(ck(7:9)-qk(1:3),fk(7:9))+...
        cross(ck(10:12)-qk(1:3),fk(10:12)))-...
        cross(qdk(1:3),diag(Ib)*qdk(1:3)));
    
    % forward euler integration of dynamics
    opti.subject_to(q(1:3,k+1)  - q(1:3,k)  == qdk(4:6) * dt(k));
    opti.subject_to(q(4:6,k+1)  - q(4:6,k)  == Binv(rpyk)*(R_body_to_world*qdk(1:3)) * dt(k));
    opti.subject_to(qdot(4:6,k+1) - qdk(4:6) == rddot * dt(k));
    opti.subject_to(qdot(1:3,k+1) - qdk(1:3) == omegaDot * dt(k));
    
    % non-negative GRF
    opti.subject_to(f_max*ones(4,1) >= fk([3 6 9 12]) >= zeros(4,1));

    % contact constraints
    R_yaw = rpyToRotMat([0 0 rpyk(3)]);
    for leg = 1:model.NLEGS
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        opti.subject_to(ck(xyz_idx(3)) >= 0);
        opti.subject_to(fk(xyz_idx(3))*ck(xyz_idx(3)) <= .01);
        if (k+1 < N)
            % no-slip constraint
            opti.subject_to(fk(xyz_idx(3))*(c(xyz_idx,k+1)-ck(xyz_idx)) <= 0.01);
            opti.subject_to(fk(xyz_idx(3))*(c(xyz_idx,k+1)-ck(xyz_idx)) >= -0.01);
        end
        
        r_hip = qk(1:3) + R_body_to_world*params.hipSrbmLocation(leg,:)';
        p_rel = (ck(xyz_idx) - r_hip);
        kin_box_x = 0.15;
        kin_box_y = 0.15;
        kin_box_z = 0.30;
        
        opti.subject_to(-kin_box_x <= p_rel(1) <= kin_box_x);
        opti.subject_to(-kin_box_y <= p_rel(2) <= kin_box_y);
        opti.subject_to(-kin_box_z <= p_rel(3) + 0.05 <= 0);
        opti.subject_to(dot(p_rel, p_rel) <= l_leg_max^2);
    end

    % friction Constraints, Eq (7k)
    opti.subject_to(fk([1 4 7 10]) <= 0.71*mu*fk([3 6 9 12]));
    opti.subject_to(fk([1 4 7 10]) >= -0.71*mu*fk([3 6 9 12]));
    opti.subject_to(fk([2 5 8 11]) <= 0.71*mu*fk([3 6 9 12]));
    opti.subject_to(fk([2 5 8 11]) >= -0.71*mu*fk([3 6 9 12]));
    
    % state & velocity bounds, Eq (7k)
    opti.subject_to(qk <= q_max);
    opti.subject_to(qk >= q_min);
    opti.subject_to(qdk <= qd_max);
    opti.subject_to(qdk >= qd_min);
    
end
%% reference trajectories
q_init_val = [0 0 0.6 0 -pi/3.2 0]';
qd_init_val = [0 0 0 1.5 -2   -3.]';

q_min_val = [-10 -10 0.075 -10 -10 -10];
q_max_val = [10 10 1.0 10 10 10];
qd_min_val = [-10 -10 -10 -40 -40 -40];
qd_max_val = [10 10 10 40 40 40];

q_term_min_val = [-10 -10 0.15 -0.1 -0.1 -10];
q_term_max_val = [10 10 5 0.1 0.1 10];
qd_term_min_val = [-10 -10 -10 -40 -40 -40];
qd_term_max_val = [10 10 10 40 40 40];

q_term_ref = [0 0 0.275, 0 0 0]';
qd_term_ref = [0 0 0, 0 0 0]';

c_init_val = repmat(q_init_val(1:3),4,1)+...
    diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -q_init_val(3)],1,4)';

c_ref = diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -0.2],1,4)';
f_ref = zeros(12,1);

QX_val = [0 0 0, 10 10 0, 10 10 10, 10 10 10]';
QX_val = zeros(12, 1);
QN_val = [0 0 100, 100 100 0, 10 10 10, 10 10 10]';
Qc_val = [0 0 0]';
Qf_val = [.001/200 .001/200 .0001/200]';

mu_val = 1;
l_leg_max_val = .35;
f_max_val = 200;

%% set parameter values
for i = 1:6
    Xref_val(i,:)   = linspace(q_init_val(i),q_term_ref(i),N);
    Xref_val(6+i,:) = linspace(qd_init_val(i),qd_term_ref(i),N);
end
for leg = 1:4
    for xyz = 1:3
        Uref_val(3*(leg-1)+xyz,:)    = Xref_val(xyz,1:end-1) + c_ref(3*(leg-1)+xyz);
        Uref_val(12+3*(leg-1)+xyz,:) = f_ref(xyz).*ones(1,N-1);
    end
end
opti.set_value(Xref, Xref_val);
opti.set_value(Uref, Uref_val);
opti.set_value(dt, dt_val);
opti.set_value(q_min, q_min_val);opti.set_value(q_max, q_max_val);
opti.set_value(qd_min, qd_min_val);opti.set_value(qd_max, qd_max_val);
opti.set_value(q_init, q_init_val);
opti.set_value(qd_init, qd_init_val);
opti.set_value(c_init, c_init_val);
opti.set_value(q_term_min, q_term_min_val);opti.set_value(q_term_max, q_term_max_val);
opti.set_value(qd_term_min, qd_term_min_val);opti.set_value(qd_term_max, qd_term_max_val);
opti.set_value(QX, QX_val);opti.set_value(QN, QN_val);
opti.set_value(Qc, Qc_val);opti.set_value(Qf, Qf_val);
opti.set_value(mu, mu_val);
opti.set_value(l_leg_max, l_leg_max_val);
opti.set_value(f_max,f_max_val);
opti.set_value(mass,mass_val);
opti.set_value(Ib,diag(Ibody_val(1:3,1:3)));
opti.set_value(Ib_inv,diag(Ibody_inv_val(1:3,1:3)));

%% initial guess
load('prevSoln.mat'); 
U_star_guess = U_star; X_star_guess = X_star; lam_g_star_guess = lam_g_star;
% opti.set_initial([U(:)],[U_star_guess(:)]);
% opti.set_initial([X(:)],[X_star_guess(:)]);
% opti.set_initial(opti.lam_g, lam_g_star);
opti.set_initial([U(:)],[Uref_val(:)]);
% opti.set_initial([X(:)],[Xref_val(:)]);

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
    'linear_solver','ma57',... % {'ma27','mumps','ma57','ma77','ma86'} % ma57 seems to work well
    'linear_system_scaling','slack-based',... {'mc19','none','slack-based'}; % Slack-based
    'linear_scaling_on_demand','yes',... % {'yes','no'};
    'max_refinement_steps',10,... % (10)
    'min_refinement_steps',1,... % (1)
    'warm_start_init_point', 'yes'); % (no)

s_opts.file_print_level = 3;
s_opts.print_level = 3;
s_opts.print_frequency_iter = 100;
s_opts.print_timing_statistics ='no';
opti.solver('ipopt',p_opts,s_opts);

 s_opts = struct('linsolver',4,... %4 works well
            'outlev', 7,...
            'strat_warm_start',0,...
            'algorithm',0,...
            'bar_murule',2,... % 5 works well
            'feastol',1e-4,...
            'tuner',0,...
            'bar_feasible',0,... %0 works well
            'bar_directinterval',10,...
            'maxit',800);%,...

% opti.solver('knitro', p_opts, s_opts);

%% solve

disp_box('Solving with Opti Stack');
tic
sol = opti.solve_limited();
toc

%% partition solution
X_star = sol.value(X);
U_star = sol.value(U);
q_star(1:6,:) = sol.value(q);
qd_star = sol.value(qdot);
lam_g_star = sol.value(opti.lam_g);
save('prevSoln.mat','X_star','U_star', 'lam_g_star');

q_foot_guess = repmat([0 -0.7 1.45]', 4, 1);

% inverse kinematics, if called
if run_IK
    for i = 1:N-1
        [x, fval, exitflag] = inverse_kinematics(U_star(1:12,i), model, q_star(1:6,i), q_foot_guess);
        if exitflag <= 0
            q_star(7:18,i) = q_foot_guess;
        end
        q_star(7:18,i) = x;
    end
    q_star(7:18,N) = q_star(7:18,N-1);
else
    q_star(7:18,:) = repmat(repmat(q_leg_home', 4, 1),1,N);
end
   
t_star = zeros(1,N);
for k = 2:N
    t_star(k) = t_star(k-1) + dt_val(1,k-1);
end
    
if(show_animation)
    showmotion(model,t_star,q_star)
end


%% plots
f_star = U_star(13:24, :); p_star = U_star(1:12, :);

if make_plots
    
    % GRFs
    figure; hold on;
    for leg = 1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        plot(t_star(1:end-1), f_star(xyz_idx(3), :));
    end
    xlabel('Time (s)'); ylabel('Force (N)');
    title('Vertical ground reaction forces');
    legend('FR', 'FL', 'BR', 'BL')
    hold off;
    
    % CoM posn
    figure; hold on;
    plot(t_star, q_star(1,:))
    plot(t_star, q_star(2,:))
    plot(t_star, q_star(3,:))
    xlabel('Time (s)'); ylabel('Position (m)');
    legend('X','Y','Z')
    title('CoM Position')
    hold off;
    
    % CoM posn
    figure; hold on;
    plot(t_star, qd_star(4,:))
    plot(t_star, qd_star(5,:))
    plot(t_star, qd_star(6,:))
    xlabel('Time (s)'); ylabel('Velocity (m/s)');
    legend('X','Y','Z')
    title('CoM Velocity')
    hold off;
    
    % orientation
    figure; hold on;
    plot(t_star, rad2deg(q_star(4,:)))
    plot(t_star, rad2deg(q_star(5,:)))
    plot(t_star, rad2deg(q_star(6,:)))
    xlabel('Time (s)'); ylabel('Orientation (degrees)');
    legend('Roll', 'Pitch', 'Yaw')
    title('CoM Orientation')
    hold off;
    
    
end
    
    
    
    
    
    
    
    
    