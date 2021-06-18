%% metadata
% Description: Piecewise polynomial evaluation test
% Author: Se Hwan Jeon

%% cleanup
clear all; clc

%% add library paths
% addpath(genpath('../../utilities/casadi/casadi_windows')); % may need to specify os directory
% addpath(genpath('../../utilities/casadi/casadi_linux'));
addpath(genpath('../../utilities'));
import casadi.*

%% test

x = casadi.MX([1,3,7,8]);

% Between 1 and 3: 0 + v + v**2
% Between 3 and 7: 3 - v - v**2
% Between 7 and 8: 4*v**2
y = casadi.MX([0 3 0; 1 -1 0; 1 -1 4]);

n = 2;

v = casadi.MX.sym('v', 1, 1);

L = low(x,v)+1;
coeff = y(:,L);

res = dot(coeff, [v.^[0:2]]');

f = casadi.Function('f',{v},{res});

% 
% for ve in [0,1,2,3,4,7,7.5,8,8.5]:
%   print(ve, f(ve), [0+ve+ve**2,3-ve-ve**2,4*ve**2] )