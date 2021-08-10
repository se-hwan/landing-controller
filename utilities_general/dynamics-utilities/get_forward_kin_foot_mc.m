function p_f = get_forward_kin_foot_mc(model, params, q)

p_f = cell(4, 1);
sideSign = [-1, 1, -1, 1];
foot_sign_convention = [1 -1 1, 1 1 1, -1 -1 1, -1 1 1];

l_1 = params.hipLocation(2);
l_2 = params.kneeLocation(3)*-1;
l_3 = params.footLocation(3)*-1;
l_4 = 0.00;

fb_state = q(1:6);
jpos = q(7:18);

R_body_to_world = rpyToRotMatTest(fb_state(4:6));

com_world = fb_state(1:3);

for leg = 1:4
    xyz_idx = 3*leg-2 : 3*leg;
    
    q_leg = jpos(xyz_idx);
    
    s1 = sin(q_leg(1)); s2 = sin(q_leg(2)); s3 = sin(q_leg(3));
    c1 = cos(q_leg(1)); c2 = cos(q_leg(2)); c3 = cos(q_leg(3));
    c23 = c2*c3 - s2*s3; s23 = s2*c3 + c2*s3;
    
    abad_world = com_world + R_body_to_world*(params.abadLocation.*foot_sign_convention(xyz_idx))';
    foot_relHip = [l_2*s2 + l_3*s23;
                   (l_1+l_4)*c1*sideSign(leg) + s1*(l_2*c2 + l_3*c23);
                   (l_1+l_4)*s1*sideSign(leg) - c1*(l_2*c2 + l_3*c23)];
    p_f{leg, 1} = abad_world + R_body_to_world*foot_relHip;
end
end

