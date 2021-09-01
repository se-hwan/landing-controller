function pf0 = get_forward_kin_knee( model, q)
% compute knee positions

%% Forward Kinematics
R_world_to_body = rpyToRotMat(q(4:6))';
for i = 1:5
    Xup{i} = zeros(6,6);
end

Xup{6} = [R_world_to_body zeros(3,3);...
    -R_world_to_body*skew(q(1:3)) R_world_to_body];

for i = 7:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
end

X0{6} = Xup{6};


for i = 7:model.NB
    X0{i} = Xup{i} * X0{model.parent(i)};               % propagate xform from origin
    [R, p0{i}] = plux_2(X0{i});                                % rotation from origin, translation from origin
    R0{i} = R';                                                % rotation **TO** origin
    Xr0{i} = [R0{i} zeros(3,3); zeros(3,3) R0{i}];             % rotation **TO** origin as spatial motion xform
end
                   
for i = 1:model.NLEGS      
    %disp(i);% loop through feet  
    hip_idx = [8 11 14 17];
    knee_idx = [9 12 15 18];
    j = hip_idx(i);                                       % body containing foot
    k = knee_idx(i);
    [~,pf0{i}] = plux_2(model.Xtree{k} * X0{j});               % origin to foot translation, world coordinates
end
