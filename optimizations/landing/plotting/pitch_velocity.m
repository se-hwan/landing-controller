clear all; clc; close all;

%% load data
addpath(genpath('../data'));
data_p30_vZ = load('pitch_30_vZ.mat');
data_p45_vX = load('pitch_45_vX.mat');
data_p45_vY = load('pitch_45_vY.mat');

data_p30_vZ = data_p30_vZ.opt_sol;
data_p45_vX = data_p45_vX.opt_sol;
data_p45_vY = data_p45_vY.opt_sol;

%% Plot Optimization results
% if not swept, velocity is -3 m/s starting from 0.6 m height

vX_plot_data = cell(4, 1);
for i = 1:4
    vX_plot_data{i} = zeros(3, numel(data_p45_vX));
end

for i = 1:numel(data_p45_vX)
    td_idx = zeros(4, 1);
    for leg = 1:4
        xyz_idx = 3*leg-2:3*leg;
        
        td_idx(leg) = find(data_p45_vX{i}.f_star(3*leg, :) > 1, 1);
        state_i = data_p45_vX{i}.X_star(:, td_idx(leg));
        p_rel_i = (data_p45_vX{i}.p_star(xyz_idx, td_idx(leg)) - data_p45_vX{i}.q_star(1:3, td_idx(leg))); 
        
        pitch_tdAngle_i = pi/4;                                   % pitch of robot at touchdown of leg i
        vel_tdAngle_i = atan(state_i(10)/state_i(end));               % velocity angle of robot at touchdown of leg i
        foot_tdAngle_i = atan(p_rel_i(1)/p_rel_i(3));
        vX_plot_data{leg}(:, i) = [pitch_tdAngle_i; vel_tdAngle_i; foot_tdAngle_i];
    end
end


figure;
hold on;
for i = 1:numel(data_p45_vX)
    plot3(vX_plot_data{1}(1, i), vX_plot_data{1}(2, i), vX_plot_data{1}(3, i), 'bo')
end
xlabel('Pitch (rad)'); ylabel('Velocity angle (rad)'); zlabel('Foot angle')
hold off;