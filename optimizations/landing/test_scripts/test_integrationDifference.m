%% metadata


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

%% tests
x_0 = 0;
v_0 = -7; 
g = 9.81;

dt = 0.03;
m = 8.5;
F = 200; 

t = 0:dt:.15;

for i = 1:length(t)
    % analytical
    x_a(i) = .5*(-g + F/m)*t(i)^2 + v_0*t(i) + x_0;

    % forward euler
    x_fe(i) = x_0 + t(i)*v_0;

    % semi-implicit euler
    v_se = v_0 + t(i)*(-g+F/m);
    x_se(i) = x_0 + t(i)*v_se;
end

%% plots
figure; hold on; 
plot(t, [x_a], 'ro-')
plot(t, [x_fe], 'go-')
plot(t, [x_se], 'bo-')
xlabel('Time')
ylabel('Position')
legend('Analytical', 'Euler', 'Semi-implicit Euler')
hold off















