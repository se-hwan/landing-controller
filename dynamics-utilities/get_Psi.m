function Psi = get_Psi(q)

% rotation matrix from RPY
rpy = q(4:6,1);
R_body_to_world = rpyToRotMat(rpy);
    
% [omega (body);vBody (body)] = phi * qdot = phi * [omega (world); vBody (world)];
Phi = [R_body_to_world' zeros(3,3);zeros(3,3) R_body_to_world'];

Psi = inv(Phi);