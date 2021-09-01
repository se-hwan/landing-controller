%% metadata
% Description:  Trajectory optimization for quadrupedal landing with single rigid body model
%               Uses Michael Posa's contact complementarity constraints with no-slip contacts
% Author:       Se Hwan Jeon

%% cleanup
rmpath(genpath('.')); % clear all previously added paths
clearvars -except training_data; clc; close all;

%% flags
show_animation = true;
run_IK = false;
make_plots = false;
make_casadi_function = true;

%% add library paths
addpath(genpath('../../../utilities_general'));
addpath(genpath('../utilities_landing'));
addpath(genpath('../codegen_casadi'));
import casadi.*

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModel(params, model);

%% timestep parameters
N = 21; 
dt_val = [0.05 repmat(0.02, 1, 15) [0.05 0.05 0.1 0.2]];

%% optimization

% q       (6) [x,y,z,roll,pitch,yaw] (floating base state)
% qdot    (6) [omega (body frame), vBody (world frame)]
% c       (NLEGS*3) (absolute world frame)
% f_grf   (NLEGS*3) (absolute world frame)

% Optimization variables
% -----------------------------------------------------------------------%
opti = casadi.Opti();
X = opti.variable(12, N);               % floating base
q       = X(1:6,:);
qdot    = X(7:12,:);
U = opti.variable(6*model.NLEGS, N-1);  % foot posns + GRFs
jpos    = opti.variable(12, N-1);         % joint positions
c     = U(1:12,:);
f_grf = U(13:24,:);

% Optimization Parameters
% -----------------------------------------------------------------------%
Xref = opti.parameter(12, N);               % floating base reference
Uref = opti.parameter(24, N-1);             % foot posn + GRF reference

dt = opti.parameter(1, N-1);                % timesteps

q_init = opti.parameter(6,1);               % initial pos. and orientation
qd_init = opti.parameter(6,1);              % initial velocities
c_init = opti.parameter(3*model.NLEGS,1);   % initial foot position

jpos_min = opti.parameter(12, 1);           % joint position limits
jpos_max = opti.parameter(12, 1);

q_term_min = opti.parameter(6,1);           % state bounds
q_term_max = opti.parameter(6,1);
qd_term_min = opti.parameter(6,1);
qd_term_max = opti.parameter(6,1);

q_min = opti.parameter(6,1);                % terminal state bounds
q_max = opti.parameter(6,1);
qd_min = opti.parameter(6,1); 
qd_max = opti.parameter(6,1);

QN = opti.parameter(12,1);                  % weighting matrices   

mu = opti.parameter();                      % friction coefficient
l_leg_max = opti.parameter();               % maximum leg length
mass = opti.parameter();                    % robot mass
Ib = opti.parameter(3,1);                   % robot inertia
Ib_inv = opti.parameter(3,1);               % robot inertia inverse

kin_box = opti.parameter(2, 1);             % foot kin. box limits (x, y)
     
%% cost function
cost = casadi.MX(0);                        % initialize cost
X_err = X(:,end)-Xref(:,end);
cost = cost + X_err'*diag(QN)*X_err;        % quadratic error terminal cost
opti.minimize(cost);                        % set objective

%% initial state constraint
opti.subject_to(q(1:6,1) == q_init);        % initial pos. + ori.
opti.subject_to(qdot(1:6,1) == qd_init);    % initial ang. vel. + lin. vel.
opti.subject_to(c(:, 1) == c_init);         % initial foot pos.

%% terminal state constraints
opti.subject_to(q(:,N) >= q_term_min);      % terminal state bounds
opti.subject_to(q(:,N) <= q_term_max);
opti.subject_to(qdot(:,N) >= qd_term_min);
opti.subject_to(qdot(:,N) <= qd_term_max);

%% general constraints

tau_leg = casadi.MX(12, N-1);

for k = 1:N-1
    % values at timestep 'k'
    qk = q(:,k);
    qdk = qdot(:,k);
    rpyk = q(4:6,k);
    ck = c(:,k);
    fk = f_grf(:,k);
    jposk = jpos(:, k);
    
    % rotation matrices
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
    opti.subject_to(qdot(4:6,k+1) - qdk(4:6) == rddot * dt(k));
    opti.subject_to(qdot(1:3,k+1) - qdk(1:3) == omegaDot * dt(k));
    opti.subject_to(q(1:3,k+1)  - q(1:3,k)  == qdk(4:6) * dt(k));
    opti.subject_to(q(4:6,k+1)  - q(4:6,k)  == Binv(rpyk)*(R_body_to_world*qdk(1:3)) * dt(k));

    % non-negative GRF
    opti.subject_to( fk([3 6 9 12]) >= zeros(4,1));

    % contact constraints
    J_f = get_foot_jacobians_mc( model, params, jposk);
    
    for leg = 1:model.NLEGS
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        opti.subject_to(ck(xyz_idx(3)) >= 0);                       % positive foot height
        opti.subject_to(fk(xyz_idx(3))*ck(xyz_idx(3)) <= .001);     % LCP constraint
        if (k+1 < N)
            % no-slip constraint
            opti.subject_to(fk(xyz_idx(3))*(c(xyz_idx,k+1)-ck(xyz_idx)) <= 0.001);
            opti.subject_to(fk(xyz_idx(3))*(c(xyz_idx,k+1)-ck(xyz_idx)) >= -0.001);
        end
        
        r_hip = qk(1:3) + R_body_to_world*params.hipSrbmLocation(leg,:)';
        p_rel = (ck(xyz_idx) - r_hip);

        kin_box_x = 0.125 + kin_box(1);
        kin_box_y = 0.125 + kin_box(2);
        kin_box_z_upper = -0.075;
        kin_box_z_lower = -0.4;
        
        sideSign = [-1, 1, -1, 1];
        
        % kinematic box and leg length constraints        
        opti.subject_to(-kin_box_x <= p_rel(1) <= kin_box_x);
        if (leg == 1 || leg == 3)
            opti.subject_to(-.05*sideSign(leg) >= p_rel(2) >= -kin_box_y);
        else
            opti.subject_to(-.05*sideSign(leg) <= p_rel(2) <= kin_box_y);
        end
        opti.subject_to(kin_box_z_lower <= p_rel(3) <= kin_box_z_upper);
        opti.subject_to(dot(p_rel, p_rel) <= l_leg_max^2);
        
        % torque constraints (J^T*F <= tau_max)
        tau_leg(xyz_idx, k) = J_f{leg}'*(-R_world_to_body*fk(xyz_idx));

        opti.subject_to(-model.tauMax(1) <= tau_leg(xyz_idx(1), k) <= model.tauMax(1));
        opti.subject_to(-model.tauMax(2) <= tau_leg(xyz_idx(2), k) <= model.tauMax(2));
        opti.subject_to(-model.tauMax(3) <= tau_leg(xyz_idx(3), k) <= model.tauMax(3));
    end
    
    % friction constraints
    opti.subject_to(fk([1 4 7 10]) <= 0.71*mu*fk([3 6 9 12]));
    opti.subject_to(fk([1 4 7 10]) >= -0.71*mu*fk([3 6 9 12]));
    opti.subject_to(fk([2 5 8 11]) <= 0.71*mu*fk([3 6 9 12]));
    opti.subject_to(fk([2 5 8 11]) >= -0.71*mu*fk([3 6 9 12]));
    
    % state & velocity bounds
    opti.subject_to(qk(3) >= q_min(3)); % only enforce z bound
    
    % joint limits and forward kinematics constraints
    pFootk = get_forward_kin_foot(model, [qk; jposk]);
    footPosk = [pFootk{1}; pFootk{2}; pFootk{3}; pFootk{4}];
    opti.subject_to(ck - footPosk >= -0.01);
    opti.subject_to(ck - footPosk <= 0.01);
    opti.subject_to(jposk >= jpos_min);
    opti.subject_to(jposk <= jpos_max);
    
%     % motor voltage limits (ignored, because motors experiencing negative work)
%     if (k+1 < N && k > 1)
%         tau_motor_des_i = tau_leg(:,i) ./ repmat(model.gr,4,1);
%         current_des_i = tau_motor_des_i ./ (1.5*repmat(model.kt, 4, 1));
%         joint_vel_i = (jposk - jpos(:, k-1))./dt_val(1);
%         back_emf_i = joint_vel_i .* repmat(model.gr, 4, 1) .* repmat(model.kt, 4, 1) * 2.0;
%         v_des_i = current_des_i .* repmat(model.Rm, 4, 1) + back_emf_i;
%         opti.subject_to(v_des_i <= model.batteryV);
%         opti.subject_to(v_des_i >= -model.batteryV);        
%     end
end

%% reference trajectories

sideSign = [1 -1 1, 1 1 1, -1 -1 1, -1 1 1];

q_init_val = [0 0 0 (.25)*(2*rand(1)-1) (pi/3)*(2*rand(1)-1) (.25)*(2*rand(1)-1)]';     % major roll
qd_init_val = [0.5*(2*rand(1,3)-1) 1*(2*rand(1, 2)-1) -4.5*rand(1)-0.5]';

for leg = 1:4
    hip_world(:, leg) = rpyToRotMat(q_init_val(4:6))*params.hipSrbmLocation(leg, :)';
end
td_hip_z = abs(min(hip_world(3,:)));

td_nom = 0.35;
z_max_td = td_nom + td_hip_z + abs(dt_val(1)*qd_init_val(6));

q_init_val(3) = z_max_td;

q_term_min_val = [-10 -10 0.15 -0.1 -0.1 -10];
q_term_max_val = [10 10 5 0.1 0.1 10];
qd_term_min_val = [-10 -10 -10 -.5 -.5 -.5];
qd_term_max_val = [10 10 10 .5 .5 .5];

q_min_val = [-10 -10 0.075 -10 -10 -10];
q_max_val = [10 10 1.0 10 10 10];
qd_min_val = [-10 -10 -10 -40 -40 -40];
qd_max_val = [10 10 10 40 40 40];

q_term_ref = [0 0 0.25, 0 0 0]';
qd_term_ref = [0 0 0, 0 0 0]';

c_init_val = zeros(12, 1);
for leg = 1:4
    xyz_idx = 3*leg-2 : 3*leg;
    p_foot_rel = sideSign(xyz_idx)'.*[0.2 0.15 -0.3]';
    c_init_val(xyz_idx) = q_init_val(1:3) + rpyToRotMat(q_init_val(4:6))*p_foot_rel;
end

q_leg_home = [0 -1.45 2.65];
q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
[~,Ibody_val] = get_mass_matrix(model, q_home, 0);
mass_val = Ibody_val(6,6);
Ibody_inv_val = inv(Ibody_val(1:3,1:3));

jpos_min_val = repmat([-pi/3, -pi/2, 0]', 4, 1);
jpos_max_val = repmat([pi/3, pi/2, 3*pi/4]', 4, 1);

v_body = rpyToRotMat(q_init_val(4:6))'*(qd_init_val(4:6));

kin_box_val = [kin_box_limits(v_body(1), 'x'); kin_box_limits(v_body(2), 'y')];

c_ref = diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.2 -0.3],1,4)';
f_ref = zeros(12,1);

QN_val = [0 0 100, 10 10 0, 10 10 10, 10 10 10]';

mu_val = 0.75;
l_leg_max_val = .4;
f_max_val = 300;

%% set parameter values
for i = 1:6
    Xref_val(i,:)   = linspace(q_init_val(i),q_term_ref(i),N);
    Xref_val(6+i,:) = linspace(qd_init_val(i),qd_term_ref(i),N);
end
for leg = 1:4
    for xyz = 1:3
        Uref_val(12+3*(leg-1)+xyz,:) = f_ref(xyz).*ones(1,N-1);
    end
end
for i = 1:N-1
    for leg = 1:4
        xyz_idx = 3*leg-2 : 3*leg;
        Uref_val(xyz_idx, i) = Xref_val(1:3, i) + rpyToRotMat(Xref_val(4:6, i))*c_ref(xyz_idx);
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
opti.set_value(QN, QN_val);
opti.set_value(mu, mu_val);
opti.set_value(l_leg_max, l_leg_max_val);
opti.set_value(mass,mass_val);
opti.set_value(Ib,diag(Ibody_val(1:3,1:3)));
opti.set_value(Ib_inv,diag(Ibody_inv_val(1:3,1:3)));
opti.set_value(jpos_min, jpos_min_val);
opti.set_value(jpos_max, jpos_max_val);
opti.set_value(kin_box, kin_box_val);

%% initial guess
f = Function.load('../codegen_casadi/landingCtrller_IPOPT.casadi');

tic
% solve problem by calling f with numerial arguments (for verification)
disp_box('Solving Problem with Solver, c code and simple bounds');
[res.x,res.f] = f(Xref_val, Uref_val,...
    dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
    q_init_val, qd_init_val, ...
    q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
    QN_val, [Xref_val(:);Uref_val(:)],...
    mu_val, l_leg_max_val, f_max_val, mass_val,...
    diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
toc

% Decompose solution
X_tmp = zeros(12, N);
U_tmp = zeros(6*model.NLEGS, N-1);

res.x = full(res.x);
X_star_guess = reshape(res.x(1:numel(X_tmp)),size(X_tmp));
U_star_guess = reshape(res.x(numel(X_tmp)+1:numel(X_tmp)+numel(U_tmp)), size(U_tmp));
opti.set_initial([U(:)],[U_star_guess(:)]);
opti.set_initial([X(:)],[X_star_guess(:)]);

jpos_guess = repmat([0, -pi/4, pi/2]', 4*(N-1), 1);
opti.set_initial(jpos(:), jpos_guess);

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
    'bound_frac',5e-3,... % (1e-2), 5e-1 works well
    'bound_push',5e-3,... % (1e-2), 5e-1 works well
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

s_opts.file_print_level = 0;
s_opts.print_level = 0;
s_opts.print_frequency_iter = 100;
s_opts.print_timing_statistics ='no';
opti.solver('ipopt',p_opts,s_opts);

%% solve
tic
disp_box('Solving with Opti Stack');
sol = opti.solve_limited();
toc

%% Generate .casadi Function
if make_casadi_function
    disp_box('Building Solver with (or without) Simple Bounds');
    nlp_opts = p_opts;
    nlp_opts.ipopt = s_opts;
    solver = nlpsol('solver','ipopt',struct('x',opti.x,'p',opti.p,'f',opti.f,'g',opti.g),nlp_opts);
    disp('Solver without Simple Bounds generated');
    
    % Generate c code for helper functions
    use_code_gen = input('Use c-generated helper functions? ');
    if use_code_gen
        make_c_helper_functions = input('Generate c code for helper functions? (takes minutes) ');
        
        nlp_opts.expand = 0;%
        c_code_folder  = '../codegen_casadi';
        c_file_name = 'landingCtrller_IPOPT_ws';
        if make_c_helper_functions % if helper functions were never generated ||~isfile('casadi_functions_gen/nlp_full_kin_stance_auto_detect.so')
            
            disp_box('Generating c code');
            
            solver.generate_dependencies([c_file_name,'.c']);% generete helper functions
            disp('Done generating .c file');
            
            disp('Compiling .c file...');
            command = ['gcc -fPIC -shared -O3 ', c_file_name,'.c -o ',c_file_name,'.so'];
            tic;
            [status1,cmdout] = system(command,'-echo'); % compile the file % takes a few minutes
            t_comp=toc;
            fprintf('Done compiling in :%.2f s\n',t_comp)
            
            if(~isfolder(c_code_folder))
                error(['Missing Folder: ',c_code_folder])
            end
            
            % move generated files
            status2 = system(['mv ',c_file_name,'.c ', c_code_folder]);
            disp([' .c file was moved to folder ', c_code_folder]);
            status3 = system(['mv ',c_file_name,'.so ', c_code_folder]);
            disp([' .so file was moved to folder ', c_code_folder]);
            
        end
        
        solver = casadi.nlpsol('solver', 'ipopt', [c_code_folder,'/',c_file_name,'.so'],nlp_opts); % load a new solver object which takes code generated dependancies
        disp('Loaded the solver with c code and simple bounds');
    end
    
    %% workaround for passing initial guess param since opti stack does not support parameterized initial guesses yet
    X_initial_guess =  casadi.MX.sym('x0',size(opti.advanced.arg.x0,1),size(opti.advanced.arg.x0,2));
    res_sym = solver('x0',X_initial_guess,'p',opti.p,'lbg',opti.lbg,'ubg',opti.ubg);
    
    % generate casadi function
    f = casadi.Function('landingCtrller_IPOPT_ws',{Xref, Uref, dt,...
        q_min, q_max, qd_min, qd_max, q_init, qd_init, ...
        q_term_min, q_term_max, qd_term_min, qd_term_max,...
        QN, X_initial_guess, ...
        mu, l_leg_max, f_max, mass, Ib, Ib_inv},{res_sym.x,res_sym.f});
    
    tic
    % solve problem by calling f with numerial arguments (for verification)
    disp_box('Solving Problem with Solver, c code and simple bounds');
    [res.x,res.f] = f(Xref_val, Uref_val,...
        dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
        q_init_val, qd_init_val, ...
        q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
        QN_val, [Xref_val(:);Uref_val(:)], ...
        mu_val, l_leg_max_val, f_max_val, mass_val,...
        diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
    toc
    
    % Decompose solution
    res.x = full(res.x);
    X_star = reshape(res.x(1:numel(X)),size(X));
    q_star(1:6,:) = X_star(1:6,:);
    q_star(7:18,:) = repmat(q_home(7:end),1,N);
    t_star = zeros(1,N);
    for k = 2:N
        t_star(k) = t_star(k-1) + dt_val(1,k-1);
    end
    
    if(show_animation==1)
        showmotion(model,t_star,q_star)
    end
    
    % Save function
    save_casadi_function = input('Save Casadi Function? ');
    if save_casadi_function
        f.save('../codegen_casadi/landingCtrller_IPOPT_ws.casadi');
        save('landingCtrller_IPOPT_ws.mat','f')
    end
    
end

%% debugging tools

% %% solve
% disp_box('Solving with Opti Stack');
% tic
% sol = opti.solve_limited();
% toc
% 
% %% partition solution
% X_star = sol.value(X);
% U_star = sol.value(U);
% q_star(1:6,:) = sol.value(q);
% qd_star = sol.value(qdot);
% lam_g_star = sol.value(opti.lam_g);
% % save('prevSoln.mat','X_star','U_star', 'lam_g_star');
% 
% q_foot_guess = repmat([0 -0.7 1.45]', 4, 1);
% 
% % inverse kinematics, if called
% if run_IK
%     for i = 1:N-1
%         [x, fval, exitflag] = inverse_kinematics(U_star(1:12,i), model, q_star(1:6,i), q_foot_guess);
%         if exitflag <= 0
%             q_star(7:18,i) = q_foot_guess;
%         end
%         q_star(7:18,i) = x;
%     end
%     q_star(7:18,N) = q_star(7:18,N-1);
% else
%     q_star(7:18,:) = repmat(repmat(q_leg_home', 4, 1),1,N);
% end
%    
% t_star = zeros(1,N);
% for k = 2:N
%     t_star(k) = t_star(k-1) + dt_val(1,k-1);
% end
%     
% if(show_animation)
%     showmotion(model,t_star,q_star)
% end
% 
% 
% %% plots
% f_star = U_star(13:24, :); p_star = U_star(1:12, :);
% 
% if make_plots
%     
%     % GRFs
%     figure; hold on;
%     for leg = 1:4
%         xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
%         plot(t_star(1:end-1), f_star(xyz_idx(3), :));
%     end
%     xlabel('Time (s)'); ylabel('Force (N)');
%     title('Vertical ground reaction forces');
%     legend('FR', 'FL', 'BR', 'BL')
%     hold off;
%     
%     % CoM posn
%     figure; hold on;
%     plot(t_star, q_star(1,:))
%     plot(t_star, q_star(2,:))
%     plot(t_star, q_star(3,:))
%     xlabel('Time (s)'); ylabel('Position (m)');
%     legend('X','Y','Z')
%     title('CoM Position')
%     hold off;
%     
%     % CoM posn
%     figure; hold on;
%     plot(t_star, qd_star(4,:))
%     plot(t_star, qd_star(5,:))
%     plot(t_star, qd_star(6,:))
%     xlabel('Time (s)'); ylabel('Velocity (m/s)');
%     legend('X','Y','Z')
%     title('CoM Velocity')
%     hold off;
%     
%     % orientation
%     figure; hold on;
%     plot(t_star, rad2deg(q_star(4,:)))
%     plot(t_star, rad2deg(q_star(5,:)))
%     plot(t_star, rad2deg(q_star(6,:)))
%     xlabel('Time (s)'); ylabel('Orientation (degrees)');
%     legend('Roll', 'Pitch', 'Yaw')
%     title('CoM Orientation')
%     hold off;
% end
    