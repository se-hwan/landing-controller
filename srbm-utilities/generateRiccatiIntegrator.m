function [RDE_step, RDE_step_fwd] = generateRiccatiIntegrator(Avbl, Bvbl)
import casadi.*

% Problem Setup
NUM_STATES = 24;
NUM_CONTROL = 12;

% States
P = SX.sym('P',NUM_STATES,NUM_STATES);

% Reference State
xref = SX.sym('xref',NUM_STATES,1);
fgrf = SX.sym('fgrf',NUM_CONTROL,1);

% Linear Dynamics
A = Avbl(xref,fgrf);
B = Bvbl(xref,fgrf);

% Weight Matrices
Q = SX.sym('Q',NUM_STATES,NUM_STATES);
R = SX.sym('R',NUM_CONTROL,NUM_CONTROL);

% Model Equations
Pdot = A'*P + P*A - P*B*(R\(B'*P)) + Q;

% Define the problem structure
rde = struct('P',P,'xref',xref,'fref',fgrf,'Q',Q,'R',R,'Pdot',Pdot);

% Continuous-time riccati differential equation as a function
f_continuous = casadi.Function('f_continuous',...
    {rde.P, rde.xref, rde.fref, rde.Q, rde.R},...
    {rde.Pdot});

% Continuous-time riccati differential equation forward integration
f_continuous_fwd = casadi.Function('f_continuous',...
    {rde.P, rde.xref, rde.fref, rde.Q, rde.R},...
    {-rde.Pdot});

% RK4 integrator that takes a single backward step
dt = SX.sym('dt');
Pf = SX.sym('Pf',NUM_STATES,NUM_STATES);
xrefk = SX.sym('xrefk',NUM_STATES,1);
fgrfk = SX.sym('fgrfk',NUM_CONTROL,1);
Qk = SX.sym('Qk',NUM_STATES,NUM_STATES);
Rk = SX.sym('Rk',NUM_CONTROL,NUM_CONTROL);

k1 = f_continuous(Pf, xrefk, fgrfk, Qk, Rk);
k2 = f_continuous(Pf + dt/2*k1, xrefk, fgrfk, Qk, Rk);
k3 = f_continuous(Pf + dt/2*k2, xrefk, fgrfk, Qk, Rk);
k4 = f_continuous(Pf + dt*k3, xrefk, fgrfk, Qk, Rk);
%P0 = Pf + dt*(k1  + 2*k2  + 2*k3  + k4 )/6;
P0 = Pf + dt*k1;
RDE_step = casadi.Function('RK4', {Pf(:), xrefk, fgrfk, Qk(:), Rk(:), dt}, {P0(:)});
RDE_step.save('casadi_functions_gen/f_RDE_step.casadi');

% RK4 integrator that takes a single forward step
k1 = f_continuous_fwd(Pf, xrefk, fgrfk, Qk, Rk);
k2 = f_continuous_fwd(Pf + dt/2*k1, xrefk, fgrfk, Qk, Rk);
k3 = f_continuous_fwd(Pf + dt/2*k2, xrefk, fgrfk, Qk, Rk);
k4 = f_continuous_fwd(Pf + dt*k3, xrefk, fgrfk, Qk, Rk);
P0 = Pf + dt*(k1  + 2*k2  + 2*k3  + k4 )/6;
RDE_step_fwd = casadi.Function('RK4', {Pf(:), xrefk, fgrfk, Qk(:), Rk(:), dt}, {P0(:)});

