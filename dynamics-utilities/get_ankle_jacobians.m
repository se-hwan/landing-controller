function [JR,JL] = get_ankle_jacobians( model, q, is_opt)

for i = 1:model.NB
    %disp(i);% loop through bodies (down)
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );           % joint xform, joint motion subspace
    Xup{i} = XJ * model.Xtree{i};                           % xform from parent
    
    if model.parent(i) == 0                              % if joint is connected to origin:
        X0{i} = Xup{i};                                     % xform from origin is xform from parent
    else                                                 % otherwise
        X0{i} = Xup{i} * X0{model.parent(i)};               % propagate xform from origin
    end
    
    R = plux_2(X0{i});                                         % rotation from origin
    R0{i} = R';                                                % rotation **TO** origin
    Xr0{i} = [R0{i} zeros(3,3); zeros(3,3) R0{i}];             % rotation **TO** origin as spatial motion xform
    
end


for i = 1:model.NB
    if(is_opt)                                            % initialize jacobian, jdqd
        J{i} = casadi.MX(zeros(6,model.NB));
    else
        J{i} = zeros(6,model.NB);
    end
    
    Xj = eye(6);                                               % from j to i (right now j = i)
    J{i}(:,i) = S{i};                                          % diagonal of jacobian is motion subspace
    j = i;
    while model.parent(j) > 0                                  % loop through i's parents (up)
        Xj = Xj * Xup{j};                                      % propagate j to i Xform
        j = model.parent(j);                                   % next j
        J{i}(:,j) = Xj * S{j};                                 % jacobian (still in i's coordinates)
    end
    
    J0{i} = Xr0{i} * J{i};                                     % jacobian (now it world coordinates)
end

JR = J0{11};%(4:6,7:11); % Right foot
JL = J0{16};%(4:6,12:16); % Left foot

% %disp(i);% loop through feet
% for i = 1:model.NLEGS
%     j = model.b_foot(i);                                       % body containing foot
%     [~,pf0{i}] = plux_2(model.Xfoot{i} * X0{j});                 % origin to foot translation, world coordinates
%     % leg jacobian (linear force component only)
%     Jf0{i} = [zeros(3,3) eye(3)] * (Xr0{j} * model.Xfoot{i} * Xr0{j}' * J0{j}); % transform world jacobian to foot
% end