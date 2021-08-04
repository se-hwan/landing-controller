function J_f = get_foot_jacobians_mc( model, params, q)

sideSign = [-1, 1, -1, 1];

l_1 = params.hipLocation(2);
l_2 = params.kneeLocation(3)*-1;
l_3 = params.footLocation(3)*-1;
l_4 = 0.004;

J_f = cell(4, 1);

for leg = 1:4
    xyz_idx = 3*leg-2 : 3*leg;
    q_leg = q(xyz_idx);
    s1 = sin(q_leg(1)); s2 = sin(q_leg(2)); s3 = sin(q_leg(3));
    c1 = cos(q_leg(1)); c2 = cos(q_leg(2)); c3 = cos(q_leg(3));
    c23 = c2*c3 - s2*s3; s23 = s2*c3 + c2*s3;

    J = [0, l_3*c23 + l_2*c2, l_3*c23;
         l_3*c1*c23 + l_2*c1*c2 - (l_1+l_4)*s1*sideSign(leg), -l_3*s1*s23 - l_2*s1*s2, -l_3*s1*s23;
         l_3*s1*c23 + l_2*c2*s1 + (l_1+l_4)*sideSign(leg)*c1, l_3*c1*s23 + l_2*c1*s2, l_3*c1*s23];
     
    J_f{leg} = J;
end

end

