%% metadata
% Description: Trajectory optimization for quadrupedal landing with single rigid body model
% Author: Se Hwan Jeon

%% cleanup
rmpath(genpath('.')); % clear all previously added paths
clear; clc; close all;

%% flags
make_casadi_function = false;
make_vbl_functions = false;
show_animation = true;
run_IK = true;

%% add library paths
% addpath(genpath('../../utilities/casadi/casadi_windows')); % may need to specify os directory
% addpath(genpath('../../utilities/casadi/casadi_linux'));
addpath(genpath('../../utilities'));
import casadi.*

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModelMC3D(params, model, 0);

%% contact schedule parameters
N = 16; % N = 11
T = 0.5; % T = 0.22
dt_val = repmat(T/(N-1),1,N-1);
cs_val = [repmat([0 0 0 0]', 1, 2) repmat([1 1 0 0]', 1, 3) repmat([1 1 1 1]', 1, 10)];
cs_TD_val = zeros(model.NLEGS,N-1);

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

cs = opti.parameter(model.NLEGS, N-1);      % contact schedule
cs_TD = opti.parameter(model.NLEGS, N-1);   % contact schedule at touchdown (?)
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
    X_err = X(:,k) - Xref(:,k);                                         % floating base error
    pf_err = repmat(X(1:3,k),model.N_GND_CONTACTS,1) + p_hip - c(:,k);  % foot position error
    U_err = U(13:24,k) - Uref(13:24,k);                                 % GRF error
    cost = cost + (X_err'*diag(QX)*X_err+...                            % sum of quadratic error costs
        pf_err'*diag(repmat(Qc,4,1))*pf_err+...
        U_err'*diag(repmat(Qf,4,1))*U_err)*dt(k);
end

X_err = X(:,end)-Xref(:,end);    % terminal cost
cost = cost + X_err'*diag(QN)*X_err;

opti.minimize(cost);             % set objective

%% initial state constraint
opti.subject_to(q(1:6,1) == q_init);        % initial pos + ori
opti.subject_to(qdot(1:6,1) == qd_init);    % initial ang. vel. + lin. vel.
opti.subject_to(c(:,1) == c_init);          % initial foot positions

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
    csk = cs(:,k);
    
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
    opti.subject_to(fk([3 6 9 12]) >= zeros(4,1));
    
    % constrain flight GRF to zero when not in contact, Eq (8a)
    opti.subject_to(fk([3 6 9 12]) <= csk.*repmat(f_max,4,1));
    
    % contact constraints
    R_yaw = rpyToRotMat([0 0 rpyk(3)]);
    for leg = 1:model.NLEGS
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        opti.subject_to(csk(leg)*ck(3*(leg-1)+3) == 0);     % foot on ground when in contact
        if (k+1 < N)                                        % no slip
            stay_on_ground = repmat(csk(leg),3,1);
            opti.subject_to(stay_on_ground.*(c(xyz_idx,k+1)-ck(xyz_idx)) == 0);
        end
        
        % kinematic Limits - applied only at touchdown
        % do these account for leg lenghts? or is that left to the timing
        % of the contact schedule?
        p_hip = qk(1:3) + R_yaw*params.hipSrbmLocation(leg,:)';
        kin_box_dim = 0.05;
        opti.subject_to( cs(leg,k)*(ck(xyz_idx(1)) - (p_hip(1) + kin_box_dim)) <= 0);
        opti.subject_to( cs(leg,k)*(ck(xyz_idx(1)) - (p_hip(1) - kin_box_dim)) >= 0);
        opti.subject_to( cs(leg,k)*(ck(xyz_idx(2)) - (p_hip(2) + kin_box_dim)) <= 0);
        opti.subject_to( cs(leg,k)*(ck(xyz_idx(2)) - (p_hip(2) - kin_box_dim)) >= 0);
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
q_init_val = [0 0 0.4 0 pi/6 0]';
qd_init_val = [0 0 0.0 1 1 -1]';

q_min_val = [-10 -10 -0 -10 -10 -10];
q_max_val = [10 10 0.4 10 10 10];
qd_min_val = [-10 -10 -10 -40 -40 -40];
qd_max_val = [10 10 10 40 40 40];

q_term_min_val = [-10 -10 0.15 -0.1 -0.1 -10];
q_term_max_val = [10 10 5 0.1 0.1 10];
qd_term_min_val = [-10 -10 -10 -40 -40 -40];
qd_term_max_val = [10 10 10 40 40 40];

q_term_ref = [0 0 0.2, 0 0 0]';
qd_term_ref = [0 0 0, 0 0 0]';

c_init_val = repmat(q_init_val(1:3),4,1)+...
    diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -q_init_val(3)],1,4)';

c_ref = diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -0.2],1,4)';
f_ref = zeros(12,1);

QX_val = [10 10 10, 10 10 10, 10 10 10, 10 10 10]';
QN_val = [0 0 100, 10 10 100, 10 10 10, 10 10 10]';
Qc_val = [0 0 0]';
Qf_val = [0.0001 0.0001 0.001]';

mu_val = 1;
l_leg_max_val = .3;
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
opti.set_value(cs, cs_val);
opti.set_value(cs_TD, cs_TD_val);
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
opti.set_initial([X(:);U(:)],[Xref_val(:);Uref_val(:)]);

%% casadi and IPOPT options
p_opts = struct('expand',true); % this speeds up ~x10
% experimental: discrete formulation of contact state
% p_opts.discrete = [zeros(12*N, 1);                  % floating base state, continuous
%                    zeros(6*model.NLEGS*(N-1),1);    % foot position + GRFS, continuous
%                    ones(model.NLEGS*(N-1), 1)];     % contact state, discrete


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

s_opts.file_print_level = 0;
s_opts.print_level = 3;
s_opts.print_frequency_iter = 100;
s_opts.print_timing_statistics ='no';
opti.solver('bonmin',p_opts,s_opts);

%% solve
disp_box('Solving with Opti Stack');
sol = opti.solve_limited();

%% partition solution
X_star = sol.value(X);
U_star = sol.value(U);
q_star(1:6,:) = sol.value(q);
q_foot_guess = repmat([0 -0.7 1.45]', 4, 1);

% post processing of foot positions during flight for inverse kinematics
for i=1:N-1
    for leg=1:4
        if cs_val(leg,i) == 0
            R_body_to_world = rpyToRotMat(q_star(4:6,i));
            U_star(3*leg-2:3*leg, i) = q_star(1:3,i) + R_body_to_world*c_ref(3*leg-2:3*leg);
        end
    end
end

% inverse kinematics, if called
if run_IK
    for i = 1:N-1
        [x, fval, exitflag] = inverse_kinematics(U_star(1:12,i), model, q_star(1:6,i), q_foot_guess, cs_val(:,i));
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

%% Plot Optimization results
if show_animation
    showmotion(model,t_star,q_star)
end

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
        c_code_folder  = 'casadi_functions_gen/c_helper_functions';
        c_file_name = 'nlp_quad_SRBM';
        if make_c_helper_functions % if helper functions were never generated ||~isfile('casadi_functions_gen/nlp_full_kin_stance_auto_detect.so')
            
            disp_box('Generating c code');
            
            solver.generate_dependencies([c_file_name,'.c']);% generete helper functions
            disp('Done generating .c file');
            
            disp('Compiling .c file (takes ~3 minutes)');
            command = ['gcc -fPIC -shared -O3 ', c_file_name,'.c -o ',c_file_name,'.so'];
            tic;
            [status1,cmdout] = system(command,'-echo'); % compile the file % takes a few minutes
            t_comp=toc;
            fprintf('Done compiling in :%.2f s\n',t_comp)
            
            if(~isfolder(c_code_folder))
                error(['Missing Folder: ',c_code_folder])
            end
            
            %move to another directory
            status2 = system(['mv ',c_file_name,'.c ', c_code_folder]);
            disp([' .c file was moved to folder ', c_code_folder]);
            status3 = system(['mv ',c_file_name,'.so ', c_code_folder]);
            disp([' .so file was moved to folder ', c_code_folder]);
            
        end
        
        solver = casadi.nlpsol('solver', 'ipopt', ['./',c_code_folder,'/',c_file_name,'.so'],nlp_opts); % load a new solver object which takes code generated dependancies
        disp('Loaded the solver with c code and simple bounds');
    end
    
    %% workaround for passing initial guess param since opti stack does not support parameterized initial guesses yet
    X_initial_guess =  casadi.MX.sym('x0',size(opti.advanced.arg.x0,1),size(opti.advanced.arg.x0,2));
    
    res_sym = solver('x0',X_initial_guess,'p',opti.p,'lbg',opti.lbg,'ubg',opti.ubg);
    % generate casadi function
    f = casadi.Function('quadSRBM',{Xref, Uref, cs, dt,...
        q_min, q_max, qd_min, qd_max, q_init, qd_init, c_init,...
        q_term_min, q_term_max, qd_term_min, qd_term_max,...
        QX, QN, Qc, Qf, X_initial_guess, mu, l_leg_max, f_max, mass,...
        Ib, Ib_inv, cs_TD},{res_sym.x,res_sym.f});
    
    % solve problem by calling f with numerial arguments (for verification)
    disp_box('Solving Problem with Solver, c code and simple bounds');
    [res.x,res.f] = f(Xref_val, Uref_val,...
        cs_val, dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
        q_init_val, qd_init_val, c_init_val,...
        q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
        QX_val, QN_val, Qc_val, Qf_val, [Xref_val(:);Uref_val(:)],...
        mu_val, l_leg_max_val, f_max_val, mass_val,...
        diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)), cs_TD_val);
    
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
        f.save('casadi_functions_gen/f_quad_SRBM.casadi');
        save('quadSRBMfun.mat','f')
    end
    
end

if make_vbl_functions
    NUM_STATES = 24;
    NUM_CONTROL = 12;
    [Avbl, Bvbl] = generateVariationalDynamics(model);
    [RDE_step, RDE_step_fwd] = generateRiccatiIntegrator(Avbl, Bvbl);
    
    % Integration Settings
    dt_riccati = 0.022;
    N_riccati = T/dt_riccati + 1;
   
    % Weight Matrices
    F = zeros(NUM_STATES,NUM_STATES);
    Fp = 1;
    Ft = 5;
    Fw = 4;
    Fv = 3;
    F(1:12,1:12) = [Fp 0 0, 0 0 0, 0 0 0, 0 0 0;...
                    0 Fp 0, 0 0 0, 0 0 0, 0 0 0;...
                    0 0 Fp, 0 0 0, 0 0 0, 0 0 0;...
                    0 0 0, Ft 0 0, 0 0 0, 0 0 0;...
                    0 0 0, 0 Ft 0, 0 0 0, 0 0 0;...
                    0 0 0, 0 0 Ft, 0 0 0, 0 0 0;...
                    0 0 0, 0 0 0, Fw 0 0, 0 0 0;...
                    0 0 0, 0 0 0, 0 Fw 0, 0 0 0;...
                    0 0 0, 0 0 0, 0 0 Fw, 0 0 0;...
                    0 0 0, 0 0 0, 0 0 0, Fv 0 0;...
                    0 0 0, 0 0 0, 0 0 0, 0 Fv 0;...
                    0 0 0, 0 0 0, 0 0 0, 0 0 Fv];
                
    Q = zeros(NUM_STATES,NUM_STATES);
    Qp = 0.25;
    Qt = 1;
    Qw = 0.5;
    Qv = 1;
    Qpw = 0;
    Qpv = 0;
    Qtw = 0;
    Q(1:12,1:12) = [Qp 0 0, 0 0 0, Qpw -Qpw 0, Qpv 0 0;...
                    0 Qp 0, 0 0 0, Qpw -Qpw 0, 0 Qpv 0;...
                    0 0 Qp, 0 0 0, Qpw -Qpw 0, 0 0 Qpv;...
                    0 0 0, Qt 0 0, Qtw 0 0, 0 0 0;...
                    0 0 0, 0 Qt 0, 0 Qtw 0, 0 0 0;...
                    0 0 0, 0 0 Qt, 0 0 Qtw, 0 0 0;...
                    Qpw Qpw Qpw, Qtw 0 0, Qw 0 0, 0 0 0;...
                    -Qpw -Qpw -Qpw, 0 Qtw 0, 0 Qw 0, 0 0 0;...
                    0 0 0, 0 0 Qtw, 0 0 Qw, 0 0 0;...
                    Qpv 0 0, 0 0 0, 0 0 0, Qv 0 0;...
                    0 Qpv 0, 0 0 0, 0 0 0, 0 Qv 0;...
                    0 0 Qpv, 0 0 0, 0 0 0, 0 0 Qv];
    
    if (min(eig(Q)) < 0)
        error('Q must be positive semi-definite')
    end
    if (min(eig(F)) < 0)
        error('Pf must be positive semi-definite')
    end
    R1 = 90.0;
    R = diag(repmat([R1 R1 R1],1,4));
    
    % Data Logging
    diag_ent = zeros(NUM_STATES,N_riccati);
    diag_ent(:,N_riccati) = diag(F);
    
    % Integrate backward in time
    P = zeros(N_riccati*NUM_STATES*NUM_STATES,1);
    P(end-NUM_STATES*NUM_STATES+1:end,1) = reshape(F,NUM_STATES*NUM_STATES,1);
    
    for k = N_riccati:-1:2
        % Sample the reference trajectory
        t_int = (k-1)*dt_riccati;
        k_opt = 1;
        while (t_int > t_star(k_opt+1) && k_opt < N-1)
            k_opt = k_opt + 1;
        end
        k_interp = (t_star(k_opt+1) - t_int) / (t_star(k_opt+1) - t_star(k_opt));
        
        xd = k_interp.*[X_star(1:12,k_opt);U_star(1:12,k_opt)] + ...
            (1-k_interp).*[X_star(1:12,k_opt+1);U_star(1:12,k_opt)];
        ud = U_star(13:24,k_opt);
        
        P_temp = RDE_step(P(NUM_STATES*NUM_STATES*(k-1)+1:NUM_STATES*NUM_STATES*k,1),...
            xd, ud, Q(:), R(:), dt_riccati);
        P(NUM_STATES*NUM_STATES*(k-2)+1:NUM_STATES*NUM_STATES*(k-1),1) = full(P_temp);
        
        diag_ent(:,k-1) = diag(full(reshape(P_temp,NUM_STATES,NUM_STATES)));
    end
    
    % Forward Integrate
    diag_ent_fwd = zeros(NUM_STATES,N_riccati);
    diag_ent_fwd(:,1) = diag_ent(:,1);
    %diag_ent_fwd(:,1) = diag(F);
    
    % Integrate backward in time
    P_fwd = zeros(N_riccati*NUM_STATES*NUM_STATES,1);
    P_fwd(1:NUM_STATES*NUM_STATES,1) = P(1:NUM_STATES*NUM_STATES,1);
    %P_fwd(1:NUM_STATES*NUM_STATES,1) = reshape(F,NUM_STATES*NUM_STATES,1);
    
    for k = 1:1:N_riccati-1
        % Sample the reference trajectory
        t_int = (k-1)*dt_riccati;
        k_opt = 1;
        while (t_int > t_star(k_opt+1) && k_opt < N-1)
            k_opt = k_opt + 1;
        end
        k_interp = (t_star(k_opt+1) - t_int) / (t_star(k_opt+1) - t_star(k_opt));
        
        xd = k_interp.*[X_star(1:12,k_opt);U_star(1:12,k_opt)] + ...
            (1-k_interp).*[X_star(1:12,k_opt+1);U_star(1:12,k_opt)];
        ud = U_star(13:24,k_opt);
        
        P_temp = RDE_step_fwd(P_fwd(NUM_STATES*NUM_STATES*(k-1)+1:NUM_STATES*NUM_STATES*k,1),...
            xd, ud, Q(:), R(:), dt_riccati);
        P_fwd(NUM_STATES*NUM_STATES*k+1:NUM_STATES*NUM_STATES*(k+1),1) = full(P_temp);
        
        diag_ent_fwd(:,k+1) = diag(full(reshape(P_temp,NUM_STATES,NUM_STATES)));
    end
    
    disp_box('Integrated the Riccati Equation');
    
    % Plot the diagonals
    figure
    subplot(4,1,1)
    plot(1:N_riccati,diag_ent(1:3,:),'r-','Linewidth',2.4)
    hold on
    plot(1:N_riccati,diag_ent_fwd(1:3,:),'k--','Linewidth',2.4)
    subplot(4,1,2)
    plot(1:N_riccati,diag_ent(4:6,:),'g','Linewidth',2.4)
    hold on
    plot(1:N_riccati,diag_ent_fwd(4:6,:),'b--','Linewidth',2.4)
    subplot(4,1,3)
    plot(1:N_riccati,diag_ent(7:9,:),'b','Linewidth',2.4)
    hold on
    plot(1:N_riccati,diag_ent_fwd(7:9,:),'g--','Linewidth',2.4)
    subplot(4,1,4)
    plot(1:N_riccati,diag_ent(10:12,:),'k','Linewidth',2.4)
    hold on
    plot(1:N_riccati,diag_ent_fwd(10:12,:),'r--','Linewidth',2.4)
    
end



