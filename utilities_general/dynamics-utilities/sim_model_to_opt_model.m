function [ model_out ] = sim_model_to_opt_model( model_in, jtype, Xtree, I, Xrot, Irot, params )
% converts Pat's simulator model to Jared's optimization model format

model_out.NB = 18;
model_out.NLEGS = 4;
model_out.gravity = [0 0 -9.81];
model_out.parent = model_in.parent;
model_out.parent(1:6) = 0:5;         % no floatbase
model_out.jtype = jtype;
model_out.jtype(1:6) = {'Px','Py','Pz','Rx','Ry','Rz'}; % no floatbase
model_out.Xtree = Xtree;
model_out.I = I;


% add feet (x4)
Xfoot_down = plux(eye(3),[0 0 -params.l2]);
model_out.Xfoot = repmat({Xfoot_down},model_out.NLEGS,1);
model_out.b_foot = [9 12 15 18];   % lower leg links

end

