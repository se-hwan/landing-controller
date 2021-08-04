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
fb_test = [0 0 0.5 0 0 0]';
q_test = repmat([0; -deg2rad(85); deg2rad(160)], 4, 1);
q_test = repmat([0; -pi/4; pi/3], 4, 1);

showmotion(model,[0, 5],repmat([fb_test; q_test], 1, 2))

% kin_box_x = -0.0;
% kin_box_y = -0.0;
% kin_box_z = -0.3;
% 
% p_FR = [kin_box_x; kin_box_y; kin_box_z];
% 
% fwd_kin = get_forward_kin_foot(model, [fb_test; q_test]);
% fwd_kin{1} = p_FR;
% p = [fwd_kin{1}' fwd_kin{2}' fwd_kin{3}' fwd_kin{4}']';

% joint_angles = quadInverseKinematics(params, fb_test, p);
% J_FR = get_foot_jacobians(model, [fb_test; joint_angles'], 0);
[JFR, JFL, JBR, JBL] = get_foot_jacobians(model, [fb_test; q_test], 0);

l_1 = params.hipLocation(2);
l_2 = params.kneeLocation(3)*-1;
l_3 = params.footLocation(3)*-1;
l_4 = 0.004;


sideSign = [-1, 1, -1, 1];

leg = 1;
    
s1 = sin(q_test(1)); s2 = sin(q_test(2)); s3 = sin(q_test(3));
c1 = cos(q_test(1)); c2 = cos(q_test(2)); c3 = cos(q_test(3));
c23 = c2*c3 - s2*s3; s23 = s2*c3 + c2*s3;

JFR
J_manual = [0, l_3*c23 + l_2*c2, l_3*c23;
            l_3*c1*c23 + l_2*c1*c2 - (l_1+l_4)*s1*sideSign(leg), -l_3*s1*s23 - l_2*s1*s2, -l_3*s1*s23;
            l_3*s1*c23 + l_2*c2*s1 + (l_1+l_4)*sideSign(leg)*c1, l_3*c1*s23 + l_2*c1*s2, l_3*c1*s23]

        
F_max_manual = inv(J_manual')*model.tauMax(1:3)
F_max_test = inv(JFR')*model.tauMax(1:3)

% F_max = (J_FR')\model.tauMax(1:3)
