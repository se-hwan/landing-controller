%% metadata
% Description:  Trajectory optimization for quadrupedal landing with single rigid body model
%               Uses Michael Posa's contact complementarity constraints with no-slip contacts
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
fb_test = [0.5 0 0.5 pi/4 pi/4 pi/6]';
q_test = repmat([0; -pi/4; pi/2], 4, 1);

showmotion(model,[0, 5],repmat([fb_test; q_test], 1, 2))

p_foot_fwd_kin = get_forward_kin_foot(model, [fb_test; q_test]);

p_foot_fwd_kin_mc = get_forward_kin_foot_mc(model, params, [fb_test; q_test])

foot_choice = 1;

foot_sign_convention = [1 -1 1, 1 1 1, -1 -1 1, -1 1 1];
hip_ref = diag(foot_sign_convention)*repmat([0.19 0.049 0],1,4)';

p_foot_fwd_kin{foot_choice};

xyz_idx = 3*(foot_choice-1)+1:3*(foot_choice-1)+3;

p_foot_test = p_foot_fwd_kin{foot_choice}(:) - hip_ref(xyz_idx);

p_x = p_foot_test(1);
p_y = p_foot_test(2);
p_z = p_foot_test(3);

l_1 = foot_sign_convention(xyz_idx(2))*0.062;
l_2 = 0.209;
l_3 = 0.195;
l_4 = foot_sign_convention(xyz_idx(2))*0.00;

%% experimental redefinitions
% hip_ref = diag(foot_sign_convention)*repmat([0.19 0.11 0],1,4)';
% l_1 = 0;

%% theta 1
th1 = atan2(p_z, p_y) + atan2(sqrt(p_y^2 + p_z^2 - l_1^2), l_1);
th1_error = rad2deg(th1) - rad2deg(q_test(1));

%% theta 2
tmp = p_y*sin(th1) - p_z*cos(th1);
A = -2*tmp*l_2; B = -2*p_x*l_2; C = l_3^2 - tmp^2 - p_x^2 - l_2^2;
th2 = atan2(B, A) + atan2(real(sqrt(A^2 + B^2 - C^2)), C);
th2_error = rad2deg(th2) - rad2deg(q_test(2));

%% theta 3
th3 = atan2(p_x - l_2*sin(th2), tmp - l_2*cos(th2)) - th2;
th3_error = rad2deg(th3) - rad2deg(q_test(3));

%% forward kinematics
R_world_to_body = rpyToRotMatTest(fb_test(4:6))';
R_body_to_world = rpyToRotMatTest(fb_test(4:6));

R_world_to_body = rpyToRotMat(fb_test(4:6))';
R_body_to_world = rpyToRotMat(fb_test(4:6));
com_world = fb_test(1:3);
hip_world = com_world + R_body_to_world*hip_ref(xyz_idx);
foot_relHip = [l_2*sin(q_test(2)) + l_3*sin(q_test(2) + q_test(3));
                          (l_1+l_4)*cos(q_test(1)) + sin(q_test(1))*(l_2*cos(q_test(2)) + l_3*cos(q_test(2)+q_test(3)));
                          (l_1+l_4)*sin(q_test(1)) - cos(q_test(1))*(l_2*cos(q_test(2)) + l_3*cos(q_test(2)+q_test(3)))];
foot_world = hip_world + R_body_to_world*foot_relHip

                      
p_foot_fwd_kin{foot_choice}
p_foot_fwd_kin_mc{foot_choice}
%%
% test = [p_foot_fwd_kin{1}; p_foot_fwd_kin{2}; p_foot_fwd_kin{3}; p_foot_fwd_kin{4}];
% out = quadInverseKinematics(params, fb_test, test)
% q_test













