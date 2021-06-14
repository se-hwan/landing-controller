function pHipY = get_forward_kin_hip( model, q)
% compute foot positions

for i = 1:14%model.NB
    %disp(i);% loop through bodies (down)
    [ XJ, ~ ] = jcalc( model.jtype{i}, q(i) );           % joint xform, joint motion subspace                                      % spatial velocity due to joint velocity
    Xup{i} = XJ * model.Xtree{i};                           % xform from parent
    
    if model.parent(i) == 0                              % if joint is connected to origin:
        X0{i} = Xup{i};                                     % xform from origin is xform from parent
    else                                                 % otherwise
        X0{i} = Xup{i} * X0{model.parent(i)};               % propagate xform from origin
    end
    
    [~, p0{i}] = plux_2(X0{i});                                % rotation from origin, translation from origin
    
end
                   
pHipY{1} = p0{9}; % right hip
pHipY{2} = p0{14};% left hip


