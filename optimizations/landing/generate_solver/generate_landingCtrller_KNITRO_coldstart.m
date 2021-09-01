%% metadata
% Description:  Trajectory optimization for quadrupedal landing with single rigid body model
%               Uses Michael Posa's contact complementarity constraints with no-slip contacts
% Author:       Se Hwan Jeon

%% cleanup
rmpath(genpath('.')); % clear all previously added paths
clear; clc; close all;

%% flags
make_casadi_function = true;
show_animation = true;
run_IK = false;
make_plots = false;

%% add library paths
addpath(genpath('../../../utilities_general'));
addpath(genpath('../codegen_casadi'));
import casadi.*

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModel(params, model);

%% contact schedule parameters
N = 21; 
T = 0.6;
dt_val = repmat(T/(N-1),1,N-1);

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
c     = U(1:12,:);
f_grf = U(13:24,:);

% Optimization Parameters
% -----------------------------------------------------------------------%
Xref = opti.parameter(12, N);           % floating base reference
Uref = opti.parameter(24, N-1);         % foot posn + GRF reference

dt = opti.parameter(1, N-1);            % timesteps

q_min = opti.parameter(6,1);
q_max = opti.parameter(6,1);
qd_min = opti.parameter(6,1); 
qd_max = opti.parameter(6,1);

q_init = opti.parameter(6,1);
qd_init = opti.parameter(6,1);

q_term_min = opti.parameter(6,1);
q_term_max = opti.parameter(6,1);
qd_term_min = opti.parameter(6,1);
qd_term_max = opti.parameter(6,1);
QN = opti.parameter(12,1);              % weighting matrices

mu = opti.parameter();                  % robot/environment parameters
l_leg_max = opti.parameter();
f_max = opti.parameter();
mass = opti.parameter();
Ib = opti.parameter(3,1);
Ib_inv = opti.parameter(3,1);

p_hip = [0.19;-0.1;-0.2;...
         0.19;0.1;-0.2;...
         -0.19;-0.1;-0.2;...
         -0.19;0.1;-0.2];
     
%% cost function
cost = casadi.MX(0);             % initialize cost
X_err = X(:,end)-Xref(:,end);    % terminal cost
cost = cost + X_err'*diag(QN)*X_err;

opti.minimize(cost);             % set objective

%% initial state constraint
opti.subject_to(q(1:6,1) == q_init);        % initial pos + ori
opti.subject_to(qdot(1:6,1) == qd_init);    % initial ang. vel. + lin. vel.

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
        opti.subject_to(fk(xyz_idx(3))*ck(xyz_idx(3)) <= .001);
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
q_init_val = [0 0 0.6 0 pi/4 -pi/6]';
qd_init_val = [0 4 5 1.3 -2 -2.]';

q_min_val = [-10 -10 0.1 -10 -10 -10];
q_max_val = [10 10 1.0 10 10 10];
qd_min_val = [-10 -10 -10 -40 -40 -40];
qd_max_val = [10 10 10 40 40 40];

q_term_min_val = [-10 -10 0.2 -0.1 -0.1 -10];
q_term_max_val = [10 10 5 0.1 0.1 10];
qd_term_min_val = [-10 -10 -10 -40 -40 -40];
qd_term_max_val = [10 10 10 40 40 40];

q_term_ref = [0 0 0.275, 0 0 0]';
qd_term_ref = [0 0 0, 0 0 0]';

c_ref = diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -0.2],1,4)';
f_ref = zeros(12,1);

QN_val = [0 0 100, 100 100 0, 10 10 10, 10 10 10]';

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
opti.set_value(q_term_min, q_term_min_val);opti.set_value(q_term_max, q_term_max_val);
opti.set_value(qd_term_min, qd_term_min_val);opti.set_value(qd_term_max, qd_term_max_val);
opti.set_value(QN, QN_val);
opti.set_value(mu, mu_val);
opti.set_value(l_leg_max, l_leg_max_val);
opti.set_value(f_max,f_max_val);
opti.set_value(mass,mass_val);
opti.set_value(Ib,diag(Ibody_val(1:3,1:3)));
opti.set_value(Ib_inv,diag(Ibody_inv_val(1:3,1:3)));

%% initial guess
opti.set_initial([U(:)],[Uref_val(:)]);
% opti.set_initial([X(:)],[Xref_val(:)]);   % generally causes difficulties converging

%% KNITRO options
p_opts = struct('expand',true); % this speeds up ~x10
s_opts = struct('linsolver',5,... %4 works well
        'maxtime_real',10.0,...
        'maxit',1500);%,...

opti.solver('knitro', p_opts, s_opts);      
disp_box('Solving with Opti Stack');
tic
sol = opti.solve_limited();
toc

%% Generate .casadi Function
if make_casadi_function
    disp_box('Building Solver with (or without) Simple Bounds');
    nlp_opts = p_opts;
    nlp_opts.knitro = s_opts;
    nlp_opts.jit=true;

    jit_options = struct('flags', '-O1', 'verbose', true, 'compiler', 'ccache gcc');
%             options = struct('jit', true, 'compiler', 'shell', 'jit_options', jit_options);
    nlp_opts.compiler='shell';
    nlp_opts.jit_options=jit_options;
    nlp_opts.jit_temp_suffix=false;
    solver = nlpsol('solver','knitro',struct('x',opti.x,'p',opti.p,'f',opti.f,'g',opti.g),nlp_opts);
    disp('Solver without Simple Bounds generated');

    %% workaround for passing initial guess param since opti stack does not support parameterized initial guesses yet
    X_initial_guess =  casadi.MX.sym('x0',size(opti.advanced.arg.x0,1),size(opti.advanced.arg.x0,2));
    res_sym = solver('x0',X_initial_guess,'p',opti.p,'lbg',opti.lbg,'ubg',opti.ubg);
    
    % generate casadi function
    f = casadi.Function('landingCtrller_KNITRO_cs',{Xref, Uref, dt,...
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
        QN_val, [Xref_val(:);Uref_val(:)],...
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
        f.save('../codegen_casadi/landingCtrller_KNITRO_cs.casadi');
%         save('landingCtrller_KNITRO_ws.mat','f')
    end
    
end


    
