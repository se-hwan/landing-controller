clear all; clc; close all;

%% load data
addpath(genpath('../data'));
addpath(genpath('../../../utilities_general'));

data_pitch_xVel = cell(13, 1);

data_p0_vX = load('pitch_0_vX.mat');
data_p5_vX = load('pitch_5_vX.mat');
data_p10_vX = load('pitch_10_vX.mat');
data_p15_vX = load('pitch_15_vX.mat');
data_p20_vX = load('pitch_20_vX.mat');
data_p25_vX = load('pitch_25_vX.mat');
data_p30_vX = load('pitch_30_vX.mat');
data_p35_vX = load('pitch_35_vX.mat');
data_p40_vX = load('pitch_40_vX.mat');
data_p45_vX = load('pitch_45_vX.mat');
data_p50_vX = load('pitch_50_vX.mat');
data_p55_vX = load('pitch_55_vX.mat');
data_p60_vX = load('pitch_60_vX.mat');

data_pitch_xVel{1} = data_p0_vX.opt_sol;
data_pitch_xVel{2} = data_p5_vX.opt_sol;
data_pitch_xVel{3} = data_p10_vX.opt_sol;
data_pitch_xVel{4} = data_p15_vX.opt_sol;
data_pitch_xVel{5} = data_p20_vX.opt_sol;
data_pitch_xVel{6} = data_p25_vX.opt_sol;
data_pitch_xVel{7} = data_p30_vX.opt_sol;
data_pitch_xVel{8} = data_p35_vX.opt_sol;
data_pitch_xVel{9} = data_p40_vX.opt_sol;
data_pitch_xVel{10} = data_p45_vX.opt_sol;
data_pitch_xVel{11} = data_p50_vX.opt_sol;
data_pitch_xVel{12} = data_p55_vX.opt_sol;
data_pitch_xVel{13} = data_p60_vX.opt_sol;


%% Parameters
p_hip = [0.19;-0.1;0;...
         0.19;0.1;-0;...
         -0.19;-0.1;0;...
         -0.19;0.1;0]; % TODO: make this an opt parameter?

%% Plot Optimization results
% if not swept, velocity is -3 m/s starting from 0.6 m height

pitch_xVel_plot_data = cell(13, 4);
for i = 1:numel(data_pitch_xVel)
    for j = 1:4
        pitch_xVel_plot_data{i, j} = zeros(3, numel(data_pitch_xVel{i}));
    end
end

for i = 1:numel(data_pitch_xVel)
    for j = 1:numel(data_pitch_xVel{i})
        td_idx = zeros(4, 1);
        for leg = 1:4
            xyz_idx = 3*leg - 2:3*leg;
            td_idx(leg) = find(data_pitch_xVel{i}{j}.f_star(3*leg, :) > 1, 1);
            
            w_R_b = rpyToRotMat(data_pitch_xVel{i}{j}.q_star(4:6, td_idx(leg)));
            state_i = data_pitch_xVel{i}{j}.X_star(:, td_idx(leg));
            % p_rel_i = data_pitch_xVel{i}{j}.p_star(xyz_idx, td_idx(leg)) + w_R_b*p_hip(xyz_idx) - data_pitch_xVel{i}{j}.q_star(1:3, td_idx(leg));
            p_rel_i = data_pitch_xVel{i}{j}.p_star(xyz_idx, td_idx(leg)) - (data_pitch_xVel{i}{j}.q_star(1:3, td_idx(leg)) + w_R_b*p_hip(xyz_idx));
            
            pitch_tdAngle_i = rad2deg(data_pitch_xVel{i}{j}.X_star(5, 1));      % pitch of robot at touchdown of leg i
            vel_tdAngle_i = rad2deg(atan(state_i(10)/state_i(end)));            % velocity angle of robot at touchdown of leg i
            foot_tdAngle_i = rad2deg(atan(p_rel_i(1)/p_rel_i(3)));
            pitch_xVel_plot_data{i, leg}(:, j) = [pitch_tdAngle_i; vel_tdAngle_i; foot_tdAngle_i];
        end
    end
end


figure;
hold on;
grid on;
for i = 3:numel(pitch_xVel_plot_data(:, 1))
    plot3(pitch_xVel_plot_data{i, 1}(1, 1:end), pitch_xVel_plot_data{i, 1}(2, 1:end), pitch_xVel_plot_data{i, 1}(3, 1:end), 'bo-')
end
title('FR Foot')
xlabel('Pitch (deg)'); ylabel('Velocity angle (deg)'); zlabel('Foot angle (deg)')
hold off;

figure;
hold on;
grid on;
for i = 3:numel(pitch_xVel_plot_data(:, 1))
    plot3(pitch_xVel_plot_data{i, 2}(1, 1:end), pitch_xVel_plot_data{i, 2}(2, 1:end), pitch_xVel_plot_data{i, 2}(3, 1:end), 'bo-')
end
title('FL Foot')
xlabel('Pitch (deg)'); ylabel('Velocity angle (deg)'); zlabel('Foot angle (deg)')
hold off;

figure;
hold on;
grid on;
for i = 3:numel(pitch_xVel_plot_data(:, 1))
    plot3(pitch_xVel_plot_data{i, 3}(1, 1:end), pitch_xVel_plot_data{i, 3}(2, 1:end), pitch_xVel_plot_data{i, 3}(3, 1:end), 'bo-')
end
title('BR Foot')
xlabel('Pitch (deg)'); ylabel('Velocity angle (deg)'); zlabel('Foot angle (deg)')
hold off;

figure;
hold on;
grid on;
for i = 3:numel(pitch_xVel_plot_data(:, 1))
    plot3(pitch_xVel_plot_data{i, 4}(1, 1:end), pitch_xVel_plot_data{i, 4}(2, 1:end), pitch_xVel_plot_data{i, 4}(3, 1:end), 'bo-')
end
title('BL Foot')
xlabel('Pitch (deg)'); ylabel('Velocity angle (deg)'); zlabel('Foot angle (deg)')
hold off;

%% Visualization
addpath(genpath('../../../utilities_general'));
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModelMC3D(params, model, 0);
t_star = linspace(0, 0.6, 41);
show_animation = true;

q_star = data_pitch_xVel{13}{6}.q_star;            % load desired trajectory to check here

if show_animation
    showmotion(model,t_star,q_star)
end


