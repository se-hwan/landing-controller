function B_dot = BmatF_dot(rpy, rpy_dot)
% Generates the B_dot matrix needed to convert euler rates into angular
% accelerations

% Equation: omega = B_dot*ThetaDot + B*ThetaDDot, where Theta are the euler angles 

theta = rpy(2);
psi = rpy(3);

thetaDot = rpy_dot(2);
psiDot = rpy_dot(3);

B_dot = [-cos(theta)*sin(psi)*psiDot - sin(theta)*thetaDot*cos(psi) -cos(psi)*psiDot 0;...
     cos(theta)*cos(psi)*psiDot - sin(theta)*thetaDot*sin(psi) -sin(psi)*psiDot 0;...
    -cos(theta)*thetaDot 0 0];