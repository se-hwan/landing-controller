clear all; clc; close all;

%%
testPoints_z = [0.5 4 0 3 0 5.5 0];
testPoints_t = [0 0.5 0.75 1.1 1.3 1.4 1.8];

pp = spline(testPoints_t, testPoints_z);

t_eval = 0:0.05:1.8;
hold on;
plot(t_eval, ppval(pp, t_eval))
plot(t_eval, polyval(pp.coefs(1,:), t_eval), 'ro-')
plot(t_eval + 0.5, polyval(pp.coefs(2,:), t_eval), 'bo-')