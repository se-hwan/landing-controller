function [ H,C,p,pf,Jf,Jdqdf,vf ] = get_simplified_3d_dynamics( model )
q = sym('q',[model.NB 1],'real');
qd = sym('qd',[model.NB 1],'real');

zero_force = repmat({zeros(6,1)},model.NLEGS,1);

[H,C,p,~,~,pf,Jf,~,~] = all_the_dynamics(model,q,qd,zero_force,1);

D_cell = { {H,C},p,pf,Jf};
D_cell = simplify_dynamics(D_cell,0);

HC = D_cell{1};
H = HC{1};
C = HC{2};
p = D_cell{2};
pf = D_cell{3};
Jf = D_cell{4};

end

