%% metadata
% Author:           Se Hwan Jeon
% Description:      Plots touchdown feet position relative to the CoM in the body frame

clc; clear all; close all; 

addpath(genpath('../..'));
import casadi.*

%% data generation

opt_fixed = 'pitch_40';
opt_sweep = 'vX';

q_init = [0 0 0.6 0 0 0]';

sweep = [-pi/3 -pi/4 -pi/6 pi/6 pi/4 pi/3];                     % sweep over initial condition
sweep = [-1.5:0.25:1.5];
p_data = cell(numel(sweep), 1);          % feet position in body frame relative to CoM
p_hip_data = cell(numel(sweep), 1);      % feet position in body frame relative to hip
v_com_data = cell(numel(sweep), 1);      % CoM velocity in body frame
dot_v_p = cell(numel(sweep), 1);         % dot product of optimized p and v_com
opt_sol = cell(1, numel(sweep));         % store optimization solutions

p_hip = [0.19;-0.1;0;...
         0.19;0.1;-0;...
         -0.19;-0.1;0;...
         -0.19;0.1;0]; % TODO: make this an opt parameter?

for i = 1:numel(sweep)
    q_init = [0 0 0.6 0 pi/4 0]';
    qd_init = [0 0 0 sweep(i) 0 -3]';
    [X_star, q_star, f_star, p_star] = eval_SRBM_CCC(q_init, qd_init);
    td_idx = zeros(4, 1);
    td(1) = find(f_star(3,:)>1, 1); td(2) = find(f_star(6,:)>1, 1);
    td(3) = find(f_star(9,:)>1, 1); td(4) = find(f_star(12,:)>1, 1);

    opt_sol{i}.X_star = X_star;
    opt_sol{i}.q_star = q_star;
    opt_sol{i}.f_star = f_star;
    opt_sol{i}.p_star = p_star;
    opt_sol{i}.td = td;
end


%% save data
save_data = input('Save optimization results? ');
if save_data
    save(['../data/',opt_fixed,'_',opt_sweep, '.mat'], 'opt_sol');
end
%% plotting

for i = 1:numel(sweep)
    p_data{i} = zeros(3, 4);
    v_com_data{i} = zeros(3, 4);
    for leg = 1:4
        xyz_idx = 3*leg-2:3*leg;
        b_R_w = rpyToRotMat(q_star(4:6, td(leg)))';
        p_data{i}(:, leg) = b_R_w*(p_star(xyz_idx, td(leg)) - q_star(1:3, td(leg))); 
        v_com_data{i}(:, leg) = b_R_w*X_star(10:12, td(leg));
    end
end

for i = 1:numel(p_data)
    leg = 4;
    for leg = 1:4
        xyz_idx = 3*leg-2:3*leg;
        p_hip_data{i}(:, leg) = p_data{i}(:, leg) - p_hip(xyz_idx);
        v_com_hat = v_com_data{i}(:, leg)./norm(v_com_data{i}(:, leg));
        p_hat = p_hip_data{i}(:, leg)./norm(p_hip_data{i}(:, leg));
        dot_v_p{i}(1, leg) = dot(v_com_hat, p_hat);
    end
end

dot_v_p

%% animation

disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModelMC3D(params, model, 0);
t_star = linspace(0, 0.6, 41);
show_animation = true;
if show_animation
    showmotion(model,t_star,q_star)
end