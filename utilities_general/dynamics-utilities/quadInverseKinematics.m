function f = quadInverseKinematics(params, fb_state, p)

% INPUTS:       params      - STRUCT, model parameters
%               fb_state    - [6 x 1], xyz posn and rpy orientation of floating base, in world frame
%               p           - [12 x 1], absolute world positions of feet

% OUTPUTS:

import casadi.*

% q = casadi.MX(zeros(12, 1));

l_2 = -params.kneeLocation(3);
l_3 = -params.footLocation(3);

foot_sign_convention = [1 -1 1, 1 1 1, -1 -1 1, -1 1 1];
hip_rel = diag(foot_sign_convention)*repmat(params.abadLocation,1,4)';
R_world_to_body = rpyToRotMat(fb_state(4:6))';
R_body_to_world = rpyToRotMat(fb_state(4:6));

for leg = 1:4
    xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
    l_1 = foot_sign_convention(xyz_idx(2))*params.hipLocation(2);
    
    p_relHip = R_world_to_body*(p(xyz_idx) - fb_state(1:3)) - hip_rel(xyz_idx);
    
    p_x = p_relHip(1); p_y = p_relHip(2); p_z = p_relHip(3);
    
    th_1 = atan2(p_z, p_y) + atan2(sqrt(p_y^2 + p_z^2 - l_1^2), l_1);
    tmp = p_y*sin(th_1) - p_z*cos(th_1);
    A = -2*tmp*l_2; B = -2*p_x*l_2; C = l_3^2 - tmp^2 - p_x^2 - l_2^2;
    if (abs(A^2 + B^2 - C^2) <= 1e-2)
        th_2 = atan2(B, A) + atan2(real(sqrt(A^2 + B^2 - C^2)), C);
    else
        th_2 = atan2(B, A) + atan2(sqrt(A^2 + B^2 - C^2), C);
    end
    
    th_3 = atan2(p_x - l_2*sin(th_2), tmp - l_2*cos(th_2)) - th_2;
    
    q(xyz_idx) = [th_1, th_2, th_3]';
end

    f = q;

end