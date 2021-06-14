function f = quadInverseKinematics(model,x)

q = [0;0;0;0;0;0;x];
pFoot_init = get_forward_kin_foot(model, q);
f = [pFoot_init{1};pFoot_init{2};pFoot_init{3};pFoot_init{4}];
