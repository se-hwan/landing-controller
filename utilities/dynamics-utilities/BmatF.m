function B = BmatF(rpy)
% Generates the B matrix needed to convert euler rates psidot, thetadot,
% and phidot into angular velocity in world frame

% Equation: omega = B*ThetaDot, where Theta is euler angles 

psi = rpy(3);
theta = rpy(2);
B = [cos(psi)*cos(theta) -sin(psi) 0;...
    cos(theta)*sin(psi) cos(psi) 0;...
    -sin(theta) 0 1];