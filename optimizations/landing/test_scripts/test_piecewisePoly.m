%% metadata
% Description: Piecewise polynomial evaluation test
% Author: Se Hwan Jeon

%% cleanup
clear all; clc; clf

%% add library paths
% addpath(genpath('../../utilities/casadi/casadi_windows')); % may need to specify os directory
% addpath(genpath('../../utilities/casadi/casadi_linux'));
addpath(genpath('../../utilities'));
import casadi.*

%% test 1: evaluate piecewise polynomial function with casadi values

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

idx = casadi.Function('idx',{v},{L});
f = casadi.Function('f',{v},{res});

% nice this shit works

%% test 2: solve for piecewise polynomial function with symbolic casadi variables

clear all; clc;
addpath(genpath('../../utilities'));
import casadi.*

opti = casadi.Opti();

time = 0:0.5:10;
dt = 0.5;
N = length(time);
poly_order = 3;
timeDuration = opti.variable(1, 4);
coeff = opti.variable(4, 4); % third order polynomials coefficients throughout timesteps

current_time = opti.parameter(1,1);
opti.set_value(current_time, dt)




v = casadi.MX.sym('v', 1, 1);

breaks = [0 timeDuration(1) sum(timeDuration(1:2)) sum(timeDuration(1:3)) sum(timeDuration(1:4))];
breaks = casadi.MX.zeros(1, 5);
for i = 1:4
    breaks(i + 1) = sum(timeDuration(1:i));
end
L = low(breaks, v) + 1;

current_poly = coeff(:, L);

res = dot(current_poly, [v.^[poly_order:-1:0]]');
idx = casadi.Function('idx',{v},{L});
f = casadi.Function('f', {v}, {res});

t = 0;
eps = 0.0001;


for i = 1:3
    opti.subject_to(polyval(coeff(i, :)', timeDuration(i)) == polyval(coeff(i+1, :)', 0))
    opti.subject_to(polyval(getDerivCoef(coeff(i, :))', timeDuration(i)) == polyval(getDerivCoef(coeff(i+1, :))', 0))
end
dt_sym = opti.parameter(1,1); opti.set_value(dt_sym, 0.5);
for i = 1:N
    t = dt_sym*i;
    current_idx = idx(t);
    opti.subject_to(coeff(current_idx, :) <= 1*ones(1,4));
end

% for i = 1:3
%     opti.subject_to(f(sum(timeDuration(1:i)) - eps) == f(sum(timeDuration(1:i)) + eps));
% end

opti.subject_to(sum(timeDuration) == 10);

% initial condition
%opti.subject_to(f(0*current_time - eps) == -5);
opti.subject_to(polyval(coeff(1, :)', 0) == -5)

opti.subject_to(polyval(coeff(end, :)', timeDuration(end)) == 0);
%opti.subject_to(coeff(4,3) == timeDuration(3))



% final condition
%opti.subject_to(f(sum(timeDuration(:)) + eps) == -10);
for i = 1:length(timeDuration)
    opti.subject_to(timeDuration(i) >= 0.5)
end


% cost = casadi.MX(0);
% for i=1:4
%     for j = 1:4
%         cost = cost + coeff(i,j)^2;
%     end
% end
% opti.minimize(cost);

opti.solver('ipopt');

%% solve
disp_box('Solving with Opti Stack');
sol = opti.solve_limited();

c_star = sol.value(coeff)
duration_star = sol.value(timeDuration)
t_star = [0 duration_star(1) sum(duration_star(1:2)) sum(duration_star(1:3)) sum(duration_star(1:4))]

pp = mkpp(t_star, c_star);
time_eval = 0:0.01:10;
plot(time_eval, ppval(pp, time_eval))
