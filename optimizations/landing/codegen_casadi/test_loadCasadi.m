%% cleanup
clear all; clc; close all;
restoredefaultpath;

%% setup
addpath(genpath('../../../utilities_general'));
import casadi.*

%% flags
run_IK = true;
show_animation = true;

%% load function
f = Function.load('landingCtrller_IPOPT.casadi');

%% parameters
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModel(params, model);

q_leg_home = [0 -1.45 2.65];
q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
[~,Ibody_val] = get_mass_matrix(model, q_home, 0);
mass_val = Ibody_val(6,6);
Ibody_inv_val = inv(Ibody_val(1:3,1:3));

N = 21; % N = 11
T = 0.6; % T = 0.22
dt_val = repmat(T/(N-1),1,N-1);

q_init_val = [0 0 0.6 0 0 0]';
qd_init_val = [0 0 0 0 0 -5]';
% q_init_val = [0 0 0.35 0 0 -pi/6]';
% qd_init_val = [0 4 5 1.3 -2 -1.]';

q_min_val = [-10 -10 0.10 -10 -10 -10];
q_max_val = [10 10 1.0 10 10 10];
qd_min_val = [-10 -10 -10 -40 -40 -40];
qd_max_val = [10 10 10 40 40 40];

q_term_min_val = [-10 -10 0.15 -0.1 -0.1 -10];
q_term_max_val = [10 10 5 0.1 0.1 10];
qd_term_min_val = [-10 -10 -10 -40 -40 -40];
qd_term_max_val = [10 10 10 40 40 40];

q_term_ref = [0 0 0.25, 0 0 0]';
qd_term_ref = [0 0 0, 0 0 0]';

c_ref = diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -0.35],1,4)';
f_ref = zeros(12,1);

QX_val = [0 0 0, 10 10 0, 10 10 10, 10 10 10]';
QX_val = zeros(12, 1);
QN_val = [0 0 100, 100 100 0, 10 10 10, 10 10 10]';
Qc_val = [0 0 0]';
Qf_val = [0.0001 0.0001 0.001]';

mu_val = .5;
l_leg_max_val = .35;
f_max_val = 500;

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

%% solve optimization

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
X = zeros(12, N);
U = zeros(6*model.NLEGS, N-1);

res.x = full(res.x);
X_star = reshape(res.x(1:numel(X)),size(X));
q_star(1:6,:) = X_star(1:6,:);
qd_star(1:6,:) = X_star(7:12,:);
q_star(7:18,:) = repmat(q_home(7:end),1,N);
U_star = reshape(res.x(numel(X)+1:numel(X)+numel(U)), size(U));
t_star = zeros(1,N);

p_star = U_star(1:12, :);
f_star = U_star(13:24, :);

for k = 2:N
    t_star(k) = t_star(k-1) + dt_val(1,k-1);
end

save('prevSoln.mat','X_star','U_star');

%% partition solution
% inverse kinematics, if called

q_foot_guess = repmat([0 -0.7 1.45]', 4, 1);
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

%% Plot Optimization results
if show_animation
    showmotion(model,t_star,q_star)
end

J_foot = cell(4, N-1);
torque = zeros(12, N-1);
for i = 1:N-1
    R_world_to_body = rpyToRotMat(q_star(4:6, i))';
    J_f = get_foot_jacobians_mc(model, params, q_star(7:18, i));
    for leg = 1:4
        xyz_idx = 3*leg-2:3*leg;
        torque(xyz_idx, i) = J_f{leg}'*(-R_world_to_body*f_star(xyz_idx, i));
    end
end

cs = zeros(4, N); cs(:, end) = ones(4, 1);
for leg = 1:4
    cs(leg, 1:end-1) = f_star(3*leg, :) > 1;
end

save('prevSoln.mat','X_star','U_star', 'q_star', 'cs');

%% plots

make_plots = true;

if make_plots
        
    % torque limits
    figure; hold on;
    plot(t_star(1:end-1), [torque(1, :);torque(4, :);torque(7, :);torque(10, :)], 'r-')
    plot(t_star(1:end-1), model.tauMax(1)*ones(1, N-1), 'r--')
    plot(t_star(1:end-1), -model.tauMax(1)*ones(1, N-1), 'r--')
    plot(t_star(1:end-1), [torque(2, :);torque(5, :);torque(8, :);torque(11, :)], 'g-')
    plot(t_star(1:end-1), model.tauMax(2)*ones(1, N-1), 'g--')
    plot(t_star(1:end-1), -model.tauMax(2)*ones(1, N-1), 'g--')
    plot(t_star(1:end-1), [torque(3, :);torque(6, :);torque(9, :);torque(12, :)], 'b-')
    plot(t_star(1:end-1), model.tauMax(3)*ones(1, N-1), 'b--')
    plot(t_star(1:end-1), -model.tauMax(3)*ones(1, N-1), 'b--')
    xlabel('Time (s)'); ylabel('Torque (Nm)')
    title("Torque Limits")
    hold off;

    % Vertical GRFs 
    figure; hold on;
    for leg = 1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        plot(t_star(1:end-1), f_star(xyz_idx(3), :));
    end
    xlabel('Time (s)'); ylabel('Force (N)');
    title('Vertical ground reaction forces');
    legend('FR', 'FL', 'BR', 'BL')
    hold off;
    
    % X GRFs
    figure; hold on;
    for leg = 1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        plot(t_star(1:end-1), f_star(xyz_idx(1), :));
    end
    xlabel('Time (s)'); ylabel('Force (N)');
    title('X ground reaction forces');
    legend('FR', 'FL', 'BR', 'BL')
    hold off;
    
    % Y GRFs
    figure; hold on;
    for leg = 1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        plot(t_star(1:end-1), f_star(xyz_idx(2), :));
    end
    xlabel('Time (s)'); ylabel('Force (N)');
    title('Y ground reaction forces');
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

    
%     % Foot locations
%     figure; hold on;
%     for leg = 1:4
%         xyz_idx = 3*leg-2 : 3*leg;
%         plot3(p_star(xyz_idx(1), :), p_star(xyz_idx(2), :), p_star(xyz_idx(3), :))
%     end
%     xlabel('x'); ylabel('y'); zlabel('z')
%     title('Foot locations')
%     hold off;
    
end
    
    
    