function [Ag,Ig] = get_centroid_momentum_matrix(H,Ic,com,rpy)

U = [eye(6) zeros(6,length(H)-6)];

% spatial Xform from G frame to body frame
R_body_to_world = rpyToRotMat(rpy);
X_G_to_i = [R_body_to_world R_body_to_world*skew(com)';...
    zeros(3,3) R_body_to_world]';

% Centroid Momentum Matrix
Ag = X_G_to_i' * U * H;

% Composite Inertia (Centroid frame)
Ig = X_G_to_i' * Ic * X_G_to_i;