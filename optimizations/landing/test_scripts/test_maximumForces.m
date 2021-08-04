%% metadata
% Description:  
% Author:       Se Hwan Jeon

%% cleanup
rmpath(genpath('.')); % clear all previously added paths
clear; clc; close all;

%% add library paths
addpath(genpath('../../../utilities_general'));
addpath(genpath('../codegen_casadi'));
import casadi.*

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModel(params, model);

%% test configurations and parameters
fb_state = [-0.19 0.049 0 0 0 0]';                   % places FR joint exactly at origin

q_test = [0; -deg2rad(5); 3*pi/4;
          0; -pi/4; pi/3;
          0; -pi/4; pi/3;
          0; -pi/4; pi/3];

showmotion(model,[0, 1],repmat([fb_state; q_test], 1, 2))

%% sweep foot positions

jpos_min = [0, -pi/3, deg2rad(15)];
jpos_max = [0, -deg2rad(5), 3*pi/4];

N_sweep = 20;

sweep_abad = linspace(jpos_min(1), jpos_max(1), N_sweep);
sweep_hip = linspace(jpos_min(2), jpos_max(2), N_sweep);
sweep_knee = linspace(jpos_min(3), jpos_max(3), N_sweep);

p_FR = [];  % foot position relative to hip (body frame)
f_z = [];   % maximum vertical force (body frame)


for i = 1:N_sweep
    for j = 1:N_sweep
        for k = 1:N_sweep
            q_eval = repmat([sweep_abad(i), sweep_hip(j), sweep_knee(k)]', 4, 1);
            p_foot = get_forward_kin_foot(model, [fb_state; q_eval]);
            p_FR = [p_FR, p_foot{1}];
            J_FR = get_foot_jacobians_mc(model, params, q_eval);
            f_max = inv(J_FR{1}')*model.tauMax(1:3);
            f_z = [f_z, f_max(3)];
        end
    end
end

%% structure data

dt = delaunayTriangulation(p_FR(1,:)',p_FR(3,:)') ;
tri = dt.ConnectivityList ;
figure; hold on;
trisurf(tri,p_FR(1,:)',p_FR(3,:)',f_z')
xlabel('x')
ylabel('z')
zlabel('F (N)')
hold off;

%% structure data


%% plotting

figure;
hold on;
for i = 1:length(f_z)
    plot3(p_FR(1, i), p_FR(3, i), f_z(i), 'b.')
end
xlabel('x')
ylabel('z')
zlabel('F (N)')





















