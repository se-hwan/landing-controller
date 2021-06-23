function [JFR,JFL,JHR,JHL] = get_foot_jacobians( model, q, is_opt)

for i = 1:model.NB
    %disp(i);% loop through bodies (down)
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );           % joint xform, joint motion subspace
    Xup{i} = XJ * model.Xtree{i};                           % xform from parent
    
    if model.parent(i) == 0                              % if joint is connected to origin:
        X0{i} = Xup{i};                                     % xform from origin is xform from parent
    else                                                 % otherwise
        X0{i} = Xup{i} * X0{model.parent(i)};               % propagate xform from origin
    end
    
end


if(is_opt == 1)                                            % initialize jacobian, jdqd
    Xc  = casadi.MX(zeros(6,6));
elseif(is_opt == 2)  
    Xc  = sym(zeros(6,6));
else
    Xc = zeros(6,6);
end

for k = 1:model.NLEGS
    if(is_opt == 1)                                            % initialize jacobian, jdqd
        Jc{k} = casadi.MX(zeros(3,model.NB));
    elseif(is_opt == 2)  
        Jc{k} = sym(zeros(3,model.NB));
    else
        Jc{k} = zeros(3,model.NB);
    end
    
    % Toe
    i = model.b_foot(k);
    R = plux_2(X0{i});                                         % rotation from origin, 
    R0 = R';                                                   % rotation **TO** origin
    
    Xc(1:3,1:3) = R0;
    Xc(4:6,1:3) = -R0*model.Xfoot{k}(4:6,1:3);
    Xc(4:6,4:6) = R0;
    Xout = Xc(4:6,:);
    
    while (i > 6)
        Jc{k}(:,i) = Xout * S{i};
        Xout = Xout * Xup{i};
        i = model.parent(i);
    end
    Jc{k}(:,1:6) = Xout;
    

end

JFR = Jc{1}(:,7:9);
JFL = Jc{2}(:,10:12);
JHR = Jc{3}(:,13:15);
JHL = Jc{4}(:,16:18);

