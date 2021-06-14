function f = humanoidInverseKinematics(model,x)

q = [0;0;0;0;0;0;x;0;0;0;0;0;0];
[pToe0,pHeel0] = get_forward_kin_toe_heel(model, q);
f = [pToe0{1};pToe0{2};pHeel0{1};pHeel0{2}];
