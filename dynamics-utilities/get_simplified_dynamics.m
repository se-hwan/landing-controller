function [ H,C,p,pf,Jf,vf ] = get_simplified_dynamics( model )

% H(q) mass matrix
% C(q,qd) bias torques
% p(q) body locations
% pf(q) foot locations
% Jf(q) foot jacobians (not spatial)
% Jdqdf(q,qd) foot (not spatial)
% vf(q,qd) foot velocity

q = sym('q',[model.NB 1],'real');
qd = sym('qd',[model.NB 1],'real');
zero_force = repmat({zeros(6,1)},model.NLEGS,1);

[H,C,p,~,~,pf,Jf,~,vf,~,~,~] = all_the_dynamics(model,q,qd,zero_force,1);
D_cell = { {H,C},p,pf,Jf,vf};
D_cell = simplify_dynamics(D_cell,1);
HC = D_cell{1};
H = HC{1};
C = HC{2};
p = D_cell{2};
pf = D_cell{3};
Jf = D_cell{4};
vf = D_cell{5};

end

