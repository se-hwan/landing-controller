clear all; clc;

%%
u_bound = linspace(50, 1000);

m = 8.25;
x0 = 0.6;
x_min = 0.1;
v_init = linspace(0, 25);

t_compress = zeros(1, length(u_bound));

v_init = -3;

for i = 1:length(t_compress)
    t_compress(i) = (-v_init + sqrt(v_init^2 - 4*(.5*(-9.81 + u_bound(i)/m)*(x0-x_min))))/(2*.5*(-9.81 + u_bound(i)/m));
end

plot(u_bound, t_compress)