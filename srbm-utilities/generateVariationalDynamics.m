function [Avbl, Bvbl] = generateVariationalDynamics(model)
import casadi.*

H = get_mass_matrix(model,zeros(18,1),0);
Ib = H(1:3,1:3);
Ib_inv = inv(Ib);
m = H(6,6);

% Reference State
p = SX.sym('p',3,1);
rpy = SX.sym('rpy',3,1);
omega = SX.sym('omega',3,1);
v = SX.sym('v',3,1);
pf = SX.sym('pf',12,1);
xref = [p;rpy;omega;v;pf];

% Error State
delta_p = SX.sym('delta_p',3,1);
delta_eta = SX.sym('delta_eta',3,1);
delta_omega = SX.sym('delta_omega',3,1);
delta_v = SX.sym('delta_v',3,1);
delta_pf = SX.sym('delta_pf',12,1);
delta_x = [delta_p;delta_eta;delta_omega;delta_v;delta_pf];

% Control
fgrf = SX.sym('fgrf',12,1);
delta_fgrf = SX.sym('delta_fgrf',12,1);

% Variational Dynamics
delta_xDot = SX(zeros(24,1));

R = rpyToRotMat(rpy)'; % R transforms body to world
delta_xDot(1:3,1) = delta_v;
delta_xDot(4:6,1) = -skew(omega)*delta_eta+delta_omega;

t1 = skew(R'*skew(pf(1:3)-p)*fgrf(1:3) + ...
    R'*skew(pf(4:6)-p)*fgrf(4:6) + ...
    R'*skew(pf(7:9)-p)*fgrf(7:9) + ...
    R'*skew(pf(10:12)-p)*fgrf(10:12)) * delta_eta;
t2a = -(skew(fgrf(1:3))*delta_pf(1:3) + ...
    skew(fgrf(4:6))*delta_pf(4:6) + ...
    skew(fgrf(7:9))*delta_pf(7:9) + ...
    skew(fgrf(10:12))*delta_pf(10:12));
t2b = skew(fgrf(1:3)+fgrf(4:6)+fgrf(7:9)+fgrf(10:12))*delta_p;
t2c = skew(pf(1:3)-p)*delta_fgrf(1:3) + ...
    skew(pf(4:6)-p)*delta_fgrf(4:6) + ...
    skew(pf(7:9)-p)*delta_fgrf(7:9) + ...
    skew(pf(10:12)-p)*delta_fgrf(10:12);
t3 = skew(Ib*omega)*delta_omega - skew(omega)*Ib*delta_omega;
delta_xDot(7:9,1) = Ib_inv*(t1 + R'*(t2a+t2b+t2c) + t3);

delta_xDot(10:12,1) = 1/m * sum(reshape(delta_fgrf,3,4),2);

delta_xDot(13:24,1) = zeros(12,1)-0.00001*delta_x(13:24,1); % small stabilizing term

% Linear System Matrices
A = jacobian(delta_xDot,delta_x);
B = jacobian(delta_xDot,delta_fgrf);

% Create Casadi Function
Avbl = casadi.Function('Avbl', {xref,fgrf}, {A});
Bvbl = casadi.Function('Bvbl', {xref,fgrf}, {B});
