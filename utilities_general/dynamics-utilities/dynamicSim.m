function xDot = dynamicSim(x,x_ref,model,tau_motor,tau_foot,f_foot,Kd_j)

xDot = zeros(2*model.NB,1);

for i = 1:model.NB
    xDot(i,1) = x(model.NB+i);
end

% Mass matrix
H = H_sym(x(1:8));
% bias torques
C = C_sym([x;zeros(6,1)]);
% Foot Jacobian
J1ek = Jf_sym1(x(1:8));
% Remove columns of J^T corresponding to tau_roll, tau_yaw, f_y
Jfek = J1ek([2 4 6],:);
% damping torque
tau_d = 1*Kd_j * x(9:16);
% Joint PD feedback torques
Kp = 50;
Kd = 1;
tau_feedback = -Kp.*(x(4:model.NB)-x_ref(4:model.NB)) - Kd.*(x(model.NB+4:end)-x_ref(model.NB+4:end));

% qddot
xDot(model.NB+1:2*model.NB,1) = H\(-Jfek.'*[tau_foot;f_foot]-C-tau_d+[0;0;0;tau_motor+tau_feedback]);