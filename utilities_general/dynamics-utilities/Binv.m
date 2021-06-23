function B = Binv(rpy)
% Generates the inverse of B matrix where the B matrix is used to convert 
% euler rates psidot, thetadot and phidot into angular velocities in world
% coordinates
%
% The inverse B matrix (calculated here) is used to convert from anuglar
% velocity in world coordinates to euler rates
%
% Suffers from singularity at theta = +- pi/2

% Equation: ThetaDot = B^(-1)*omega

psi = rpy(3);
theta = rpy(2);
B = [cos(psi)/cos(theta) sin(psi)/cos(theta) 0;
    -sin(psi) cos(psi) 0;
    cos(psi)*tan(theta) sin(psi)*tan(theta) 1];