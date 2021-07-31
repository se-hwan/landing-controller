%% metadata
% Description:      Samples from 9D initial falling conditions and runs a
%                   trajectory optimization for each initial point. Saves
%                   the solutions to ../data
% Author:           Se Hwan Jeon

%% cleanup
restoredefaultpath;
clear; clc; 

%% setup
addpath(genpath('../../../utilities_general'));
import casadi.*

%% parameters
generate_optimization_parameters();
load('optimization_parameters.mat');

%% load solvers
f_cs = Function.load('../codegen_casadi/landingCtrller_IPOPT.casadi');
f_ws = Function.load('../codegen_casadi/landingCtrller_IPOPT_ws.casadi');

%% sample points for initial falling conditions
q_init_sample = zeros(9, 5);
q_init_sample(1, :) = linspace(-pi/3, pi/3, 5);       % roll
q_init_sample(2, :) = linspace(-pi/4, pi/3, 5);       % pitch
q_init_sample(3, :) = linspace(-pi/2, pi/2, 5);       % yaw
q_init_sample(4, :) = linspace(-3, 3, 5);             % x angular velocity
q_init_sample(5, :) = linspace(-3, 3, 5);             % y angular velocity
q_init_sample(6, :) = linspace(-3, 3, 5);             % z angular velocity
q_init_sample(7, :) = linspace(-2.5, 2.5, 5);         % x linear velocity
q_init_sample(8, :) = linspace(-2.5, 2.5, 5);         % y linear velocity
q_init_sample(9, :) = linspace(-3, -1, 5);          % z linear velocity

%% data matrix
training_data_inputs = zeros(1, 9);
training_data_outputs = zeros(1, 624);

for cntr = 1:5
    for j = 1:5
        %% set initial falling conditions
        q_init_val = [0 0 0.6 0 q_init_sample(2,j) 0]';
        qd_init_val = [0 0 0 0 0 q_init_sample(9,cntr)]';

        q_term_val = [0 0 0.25, 0 0 0]';
        qd_term_val = [0 0 0, 0 0 0]';

        %% generate reference trajectory (linearly interpolated, "cold start")
        [Xref_val, Uref_val] = generate_reference_traj(q_init_val, qd_init_val, q_term_val, qd_term_val, N);

        tic
        [res_cs.x,res_cs.f] = f_cs(Xref_val, Uref_val,...
            dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
            q_init_val, qd_init_val, ...
            q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
            QN_val, [Xref_val(:);Uref_val(:)],...
            mu_val, l_leg_max_val, f_max_val, mass_val,...
            Ibody_diag, Ibody_inv_diag);
        toc

        % Decompose solution
        X = zeros(12, N);
        U = zeros(6*model.NLEGS, N-1);

        res_cs.x = full(res_cs.x);
        X_star_cs = reshape(res_cs.x(1:numel(X)),size(X));
        U_star_cs = reshape(res_cs.x(numel(X)+1:numel(X)+numel(U)), size(U));

        tic
        [res_ws.x,res_ws.f] = f_ws(X_star_cs, U_star_cs,...
            dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
            q_init_val, qd_init_val, ...
            q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
            QN_val, [X_star_cs(:); U_star_cs(:)],...
            mu_val, l_leg_max_val, f_max_val, mass_val,...
            Ibody_diag, Ibody_inv_diag);
        toc

        disp("Falling conditions: ")
        disp("q:     " + num2str(q_init_val'))
        disp("q_dot: "+ num2str(qd_init_val')) 
        
        res_ws.x = full(res_ws.x);
        X_star = reshape(res_ws.x(1:numel(X)),size(X));
        q_star(1:6,:) = X_star(1:6,:);
        q_star(7:18,:) = repmat(q_home(7:end),1,N);
        U_star = reshape(res_ws.x(numel(X)+1:numel(X)+numel(U)), size(U));
        t_star = zeros(1,N);

        p_star = U_star(1:12, :);
        f_star = U_star(13:24, :);

        for k = 2:N
            t_star(k) = t_star(k-1) + dt_val(1,k-1);
        end

        q_foot_guess = repmat([0 -0.7 1.45]', 4, 1);
        for i = 1:N-1
            [x, fval, exitflag] = inverse_kinematics(U_star(1:12,i), model, q_star(1:6,i), q_foot_guess);
            if exitflag <= 0
                q_star(7:18,i) = q_foot_guess;
            end
            q_star(7:18,i) = x;
        end
        q_star(7:18,N) = q_star(7:18,N-1);
        showmotion(model,t_star,q_star)
        pause;
    end
end
% save([ , , '.mat'], 'X_star','U_star','t_star','q_star')







