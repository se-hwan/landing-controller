function [JtoeR,JtoeL,JheelR,JheelL] = get_toe_heel_approx_jacobians( model, task, q_in, is_opt)

if(is_opt == 1)                                            
    q  = casadi.MX(tsak.q_home);
elseif(is_opt == 2)  
    q  = sym(task.q_home);
else
    q = task.q_home;
end
q([9 10 14 15]) = q_in([9 10 14 15]); % knee positions

for i = 1:model.NB
    %disp(i);% loop through bodies (down)
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );           % joint xform, joint motion subspace
    Xup{i} = XJ * model.Xtree{i};                           % xform from parent
    
    if model.parent(i) == 0                              % if joint is connected to origin:
        X0{i} = Xup{i};                                     % xform from origin is xform from parent
    else                                                 % otherwise
        X0{i} = Xup{i} * X0{model.parent(i)};               % propagate xform from origin
    end
    
    %R = plux_2(X0{i});                                         % rotation from origin
    %R0{i} = R';                                                % rotation **TO** origin
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
        Jc_toe{k} = casadi.MX(zeros(3,model.NB));
        Jc_heel{k} = casadi.MX(zeros(3,model.NB));
    elseif(is_opt == 2)  
        Jc_toe{k} = sym(zeros(3,model.NB));
        Jc_heel{k} = sym(zeros(3,model.NB));
    else
        Jc_toe{k} = zeros(3,model.NB);
        Jc_heel{k} = zeros(3,model.NB);
    end
    
    % Toe
    i = model.b_toe(k);
    R = plux_2(X0{i});                                         % rotation from origin, 
    R0 = R';                                                   % rotation **TO** origin
    
    Xc(1:3,1:3) = R0;
    Xc(4:6,1:3) = -R0*model.Xtoe{k}(4:6,1:3);
    Xc(4:6,4:6) = R0;
    Xout = Xc(4:6,:);
    
    while (i > 6)
        Jc_toe{k}(:,i) = Xout * S{i};
        Xout = Xout * Xup{i};
        i = model.parent(i);
    end
    Jc_toe{k}(:,1:6) = Xout;
    
    % Heel
    i = model.b_heel(k);
    R = plux_2(X0{i});                                         % rotation from origin, 
    R0 = R';                                                   % rotation **TO** origin
    
    Xc(1:3,1:3) = R0;
    Xc(4:6,1:3) = -R0*model.Xheel{k}(4:6,1:3);
    Xc(4:6,4:6) = R0;
    Xout = Xc(4:6,:);
    
    while (i > 6)
        Jc_heel{k}(:,i) = Xout * S{i};
        Xout = Xout * Xup{i};
        i = model.parent(i);
    end
    Jc_heel{k}(:,1:6) = Xout;

end

JtoeR = Jc_toe{1}(:,7:11);
JtoeL = Jc_toe{2}(:,12:16);
JheelR = Jc_heel{1}(:,7:11);
JheelL = Jc_heel{2}(:,12:16);

