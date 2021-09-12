%% metadata
% Description:  Trajectory optimization for quadrupedal landing with single rigid body model
%               Uses Michael Posa's contact complementarity constraints with no-slip contacts
% Author:       Se Hwan Jeon

%% cleanup
rmpath(genpath('.')); % clear all previously added paths
clearvars -except training_data; clc; close all;

%% flags
show_animation = false;
show_plots = false;

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
dt_val = [0.05 repmat(0.02, 1, 15) [0.05 0.05 0.15 0.2]];

%% load function
f_ipopt_SRBM = Function.load('../codegen_casadi/landingCtrller_IPOPT.casadi');
f_knitro = Function.load('../codegen_casadi/landingCtrller_KNITRO.casadi');
f_knitro_ws = Function.load('../codegen_casadi/landingCtrller_KNITRO_ws.casadi');

%% main optimization loop
N_samples = 5000; % number of optimizations

for cntr = 1:N_samples

    %% reference trajectories

    sideSign = [1 -1 1, 1 1 1, -1 -1 1, -1 1 1];

    q_init_val = [0 0 0 (pi/6)*randn(1) (pi/4)*randn(1) (pi/6)*randn(1)]';     % major pitch
    qd_init_val = [0.5*randn(1,3) 1.5*(2*rand(1, 2)-1) -4*rand(1)-2]';
    
%     q_init_val = [0 0 0 0 pi/3 0]'; qd_init_val = [0 0 0 0 0 -3]';
%     q_init_val = [0 0 0 .25 pi/6 0.1]';     % major yaw 
%     qd_init_val = [0.5*(2*rand(1,3)-1) 1*(2*rand(1, 2)-1) 0]';

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

    q_term_ref = [0 0 0.225, 0 0 0]';
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
    f_max_val = 400;

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
    
    U_guess = zeros(24, 20);
    
    %% initial guess (SRBM LCP solution)
    % Decompose solution
    X_tmp = zeros(12, N);
    U_tmp = zeros(6*model.NLEGS, N-1);
    jpos_tmp = zeros(3*model.NLEGS, N-1);

    tic
    % solve problem by calling f with numerial arguments (for verification)
    disp_box('Solving coldstart...');
    [res.x,res.f] = f_ipopt_SRBM(Xref_val, Uref_val,...
        dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
        q_init_val, qd_init_val, ...
        q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
        QN_val, [Xref_val(:); Uref_val(:)],...
        mu_val, l_leg_max_val, f_max_val, mass_val,...
        diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
    toc

    res.x = full(res.x);
    X_star_guess = reshape(res.x(1:numel(X_tmp)),size(X_tmp));
    U_star_guess = reshape(res.x(numel(X_tmp)+1:numel(X_tmp)+numel(U_tmp)), size(U_tmp));
    jpos_guess = repmat([0, -pi/4, pi/2]', 4*(N-1), 1);
    
    
    %% solve
    disp_box('Solving KD problem...');
    tic
    [res.x,res.f] = f_knitro(Xref_val, Uref_val,...
            dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
            q_init_val, qd_init_val, c_init_val, ...
            q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
            QN_val, [X_star_guess(:); U_star_guess(:); jpos_guess],...
            jpos_min_val, jpos_max_val, kin_box_val, mu_val, l_leg_max_val, mass_val,...
            diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
    toc

    %% partition solution
    res.x = full(res.x);
    X_star_guess = reshape(res.x(1:numel(X_tmp)),size(X_tmp));
    U_star_guess = reshape(res.x(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
    jpos_star_guess = reshape(res.x(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(U_tmp) + numel(jpos_tmp)), size(jpos_tmp));

        res.x = full(res.x);
    X_star = reshape(res.x(1:numel(X_tmp)),size(X_tmp));
    U_star = reshape(res.x(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
    jpos_star = reshape(res.x(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(U_tmp) + numel(jpos_tmp)), size(jpos_tmp));

    q_star(1:6,:) = X_star(1:6,:); qd_star(1:6,:) = X_star(7:12,:);
    q_star(7:18,1:end-1) = jpos_star; q_star(7:18, end) = q_star(7:18, end-1);

    f_star = U_star(13:24, :); p_star = U_star(1:12, :);

    t_star = zeros(1,N);
    for k = 2:N
        t_star(k) = t_star(k-1) + dt_val(1,k-1);
    end
    if show_animation
        showmotion_floatingBase(model,t_star(1:end-1),q_star(:,1:end-1))
        plot_results(model, params, t_star, X_star, U_star, jpos_star);
    end
    
    %% resolve with previous solution

    disp_box('Final KD solve with warm start...');
    tic
    try
        [res.x,res.f] = f_knitro_ws(Xref_val, Uref_val,...
                dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
                q_init_val, qd_init_val, c_init_val, ...
                q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
                QN_val, [X_star_guess(:); U_star_guess(:); jpos_star_guess(:); ],...
                jpos_min_val, jpos_max_val, kin_box_val, mu_val, l_leg_max_val, mass_val,...
                diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
        t_ws_solve = toc;

        disp("Initial position conditions: [" + num2str(q_init_val') + "]");
        disp("Initial velocity conditions: [" + num2str(qd_init_val') + "]");

        res.x = full(res.x);
        X_star = reshape(res.x(1:numel(X_tmp)),size(X_tmp));
        U_star = reshape(res.x(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
        jpos_star = reshape(res.x(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(U_tmp) + numel(jpos_tmp)), size(jpos_tmp));

        q_star(1:6,:) = X_star(1:6,:); qd_star(1:6,:) = X_star(7:12,:);
        q_star(7:18,1:end-1) = jpos_star; q_star(7:18, end) = q_star(7:18, end-1);

        f_star = U_star(13:24, :); p_star = U_star(1:12, :);

        t_star = zeros(1,N);
        for k = 2:N
            t_star(k) = t_star(k-1) + dt_val(1,k-1);
        end

        %% visualization
        if show_animation
            showmotion_floatingBase(model,t_star(1:end-1),q_star(:,1:end-1))
        end
        if show_plots
            plot_results(model, params, t_star, X_star, U_star, jpos_star);
        end

        %% check for ground penetration and resteps
        is_valid = check_solution(model, q_star, f_star)
%     disp(is_valid);
    catch
        disp("OPTIMIZATION FAILED TO CONVERGE");
        is_valid = false;
    end
    
    %% training data generation
    disp("Iteration: " + num2str(cntr));
    if (is_valid)
        disp("Trajectory saved!. Continuing...")
        if (exist('training_data','var') == 1)
            training_data.input(:, end + 1) = [q_init_val(4:6); qd_init_val(:)];
            training_data.output(:, end + 1) = [X_star(:); U_star(:); jpos_star(:)];
        else
            if (isfile("training_data_landing.mat"))
                load("training_data_landing.mat")
                training_data.input(:, end + 1) = [q_init_val(4:6); qd_init_val(:)];
                training_data.output(:, end + 1) = [X_star(:); U_star(:); jpos_star(:)];
            else
                training_data = struct();
                training_data.input(:, 1) = [q_init_val(4:6); qd_init_val(:)];
                training_data.output(:, 1) = [X_star(:); U_star(:); jpos_star(:)];
            end
        end
    else
        disp("Trajectory not saved. Continuing...")
    end
    save("training_data_landing.mat", 'training_data', '-V7.3');

end
    
%% show distribution of training data
show_training_distribution = input("Show training data input distribution? ");
if show_training_distribution == 1
    t = tiledlayout(3, 3);
    title("Training data input distributions")
    
    nexttile; histogram(training_data.input(1,:));
    xlabel('Roll'); ylabel('Frequency');
    
    nexttile; histogram(training_data.input(2,:));
    xlabel('Pitch'); ylabel('Frequency');
    
    nexttile; histogram(training_data.input(3,:));
    xlabel('Yaw'); ylabel('Frequency');
    
    nexttile; histogram(training_data.input(4,:));
    xlabel('\omega_x'); ylabel('Frequency');
    
    nexttile; histogram(training_data.input(5,:));
    xlabel('\omega_y'); ylabel('Frequency');
    
    nexttile; histogram(training_data.input(6,:));
    xlabel('\omega_z'); ylabel('Frequency');
    
    nexttile; histogram(training_data.input(7,:));
    xlabel('V_x'); ylabel('Frequency');
    
    nexttile; histogram(training_data.input(8,:));
    xlabel('V_y'); ylabel('Frequency');
    
    nexttile; histogram(training_data.input(9,:));
    xlabel('V_z'); ylabel('Frequency');
end

%% test samples of training data

sample_idx = randi(1658);

res.x = training_data.output(:, sample_idx);
X_star = reshape(res.x(1:numel(X_tmp)),size(X_tmp));
U_star = reshape(res.x(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
jpos_star = reshape(res.x(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(U_tmp) + numel(jpos_tmp)), size(jpos_tmp));

q_star(1:6,:) = X_star(1:6,:); qd_star(1:6,:) = X_star(7:12,:);
q_star(7:18,1:end-1) = jpos_star; q_star(7:18, end) = q_star(7:18, end-1);

f_star = U_star(13:24, :); p_star = U_star(1:12, :);

t_star = zeros(1,N);
for k = 2:N
    t_star(k) = t_star(k-1) + dt_val(1,k-1);
end

showmotion_floatingBase(model,t_star(1:end-1),q_star(:,1:end-1))
plot_results(model, params, t_star, X_star, U_star, jpos_star);

