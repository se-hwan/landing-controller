%% metadata
% Description: MINLP test with Casadi
% Author: Se Hwan Jeon

%% cleanup
clc; clear all; 

%% add library paths
addpath(genpath('casadi/casadi_windows')); % choose appropriate OS
% addpath(genpath('casadi/casadi_linux'));
import casadi.*

%% problem
x = MX.sym('x');
y = MX.sym('y');

prob = struct;
prob.x = [x;y];
prob.f = (x-1.3)^2+(y-1.4)^2;

opts = struct;
opts.discrete = [true false];
opts.discrete = [0 1];

solver = nlpsol('solver', 'bonmin', prob, opts);

res = solver('x0',0);
res.x