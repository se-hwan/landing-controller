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
params.model()
model  = get_robot_model(params);
model  = buildShowMotionModel(params, model);

%% mass test
q_leg_home = [0 -1.45 2.65];
q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
[~,Ibody_val] = get_mass_matrix(model, q_home, 0);
mass_val = Ibody_val(6,6)
Ibody_inv_val = inv(Ibody_val(1:3,1:3));

%% tests
r = [1 0 0]';
O = [0 0 0]';

rotation = [pi/4 pi/4 0]; % roll, pitch, yaw 
R_zyx = rpyToRotMat(rotation);
R_xyz = rpyToRotMatTest(rotation);

r_zyx = R_zyx*r;
r_xyz = R_xyz*r;

%% visualization

baseline = [O r];
xyzRot = [O r_xyz];
zyxRot = [O r_zyx];

figure; hold on; grid on; axis equal;

plot3(xyzRot(1,:), xyzRot(2,:), xyzRot(3,:), 'ro-' )
plot3(zyxRot(1,:), zyxRot(2,:), zyxRot(3,:), 'bo-' )
plot3(baseline(1,:), baseline(2,:), baseline(3,:), 'go--' )
plot3([0 0], [0 1], [0 0], 'ko--' )
plot3([0 1], [0 0], [0 0], 'ko--' )
plot3([0 0], [0 0], [0 1], 'ko--' )
xlabel('x'); ylabel('y'); zlabel('z')
legend('xyz','zyx','original')

fb_test = [0 0 0.5 rotation]';
% fb_test = [0; 0; 0.5; -0.4377; -0.6194; -0.3749];
q_test = [0; -0.8; 1.6; zeros(9, 1)];

showmotion_floatingBase(model,[0, 1],repmat([fb_test; q_test], 1, 2))

















