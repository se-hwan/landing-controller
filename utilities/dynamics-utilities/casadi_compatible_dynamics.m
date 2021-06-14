function [ H,C,p0,J0,Jdqd,pf0,Jf0,Jdfqd0,vvf ] = casadi_compatible_dynamics( model, q, qd, f, is_symbolic )
% compute H,C,body position, body jacobians, body Jdqd's, foot positions,
% foot Jdqd's, foot velocities

% f is spatial force on bodies containing feet
% must be expressed wrt BODY FRAME!

if(is_symbolic)
    disp(['Calculating symbolic dynamics for system with ', num2str(length(q)), ' variables']);
end

a_grav = [0 0 0 0 0 -9.81]';                                % spatial acceleration of gravity


for i = 1:model.NLEGS                                       % loop through feet
    foot_body = model.b_foot(i);                            % rigid body containing foot
    f_ext{foot_body} = f{i};                                % spatial force on body due to foot force
end


if(is_symbolic)
    disp('NEWTON-EULER, first pass');
end
for i = 1:model.NB
    %disp(i);% loop through bodies (down)
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );           % joint xform, joint motion subspace
    vJ = S{i}*qd(i);                                        % spatial velocity due to joint velocity
    Xup{i} = XJ * model.Xtree{i};                           % xform from parent
    
    if model.parent(i) == 0                              % if joint is connected to origin:
        X0{i} = Xup{i};                                     % xform from origin is xform from parent
        v{i} = vJ;                                          % spatial velocity is joint velocity
        avp{i} = Xup{i} * -a_grav;                          % spatial acceleration is gravity
        a_no_qdd{i} = zeros(6,1);                           % acceleration due to joint velocity is 0
    else                                                 % otherwise
        X0{i} = Xup{i} * X0{model.parent(i)};               % propagate xform from origin
        v{i} = Xup{i}*v{model.parent(i)} + vJ;              % propogate velocity
        avp{i} = Xup{i}*avp{model.parent(i)} + crm(v{i})*vJ;   % propogate accelerations (X*vd + Xd*v)
        a_no_qdd{i} = Xup{i}*a_no_qdd{model.parent(i)} + crm(v{i})*vJ;  
    end
    
    fvp{i} = model.I{i}*avp{i} + crf(v{i})*model.I{i}*v{i};    % intermediate spatial force for bias torque
    [R, p0{i}] = plux_2(X0{i});                                  % rotation from origin, translation from origin
    R0{i} = R';                                                % rotation **TO** origin
    Xr0{i} = [R0{i} zeros(3,3); zeros(3,3) R0{i}];             % rotation **TO** origin as spatial motion xform
    
end


if(is_symbolic)
    disp('NEWTON-EULER, external forces');
end
for i = 1:model.NB    
    %disp(i);% loop through bodies
    if length(f_ext{i}) > 0    
        if(any(f_ext{i})) % huge time savings here!!!
            fvp{i} = fvp{i} - X0{i}' \ f_ext{i};  % transform it into body frame with spatial force xform
        end
    end                                                    %   remove it from fvp
end


if(is_symbolic)
    disp('NEWTON-EULER, calc_C');
else
C = casadi.MX(model.NB,1); 
end
for i = model.NB:-1:1  
    %disp(i);% loop through bodies backward (final pass of newton euler)
    C(i,1) = S{i}' * fvp{i};   
    if model.parent(i) ~= 0
        fvp{model.parent(i)} = fvp{model.parent(i)} + Xup{i}'*fvp{i};  
    end
end


if(is_symbolic)
    disp('CRBA, inertia calc');
end

IC = model.I;
for i = model.NB:-1:1 
    %disp(i);% calc composite inertia (from spatial v2 CRBA implementation)
    if model.parent(i) ~= 0
        IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}'*IC{i}*Xup{i}; 
    end
end


if is_symbolic
    H = sym(zeros(model.NB));
else
    H = casadi.MX(model.NB,model.NB);
end

if(is_symbolic)
    disp('CRBA, calc_H');
end
for i = 1:model.NB
    %disp(i);
    if(is_symbolic)                                            % initialize jacobian, jdqd
        J{i} = sym(zeros(6,model.NB));
        Jdqd{i} = sym(zeros(6,1));
    else
        J{i} = casadi.MX(6,model.NB);
        Jdqd{i} = casadi.MX(6,1);
    end
    
    Xj = eye(6);                                               % from j to i (right now j = i)
    J{i}(:,i) = S{i};                                          % diagonal of jacobian is motion subspace
    Jdqd{i} = Xr0{i} * a_no_qdd{i};                            % Jdqd is spatial acceleration due to qd, rotated to world
    v0 = Xr0{i} * v{i};                                        % spatial velocity, rotated to world
    Jdqd{i}(4:6) = Jdqd{i}(4:6) + skew(v0(1:3)) * v0(4:6);     % spatial -> normal acceleration
    
    fh = IC{i} * S{i};                                         % CRBA
    H(i,i) = S{i}' * fh;                                       % CRBA
    j = i;                  
    while model.parent(j) > 0                                  % loop through i's parents (up)
        Xj = Xj * Xup{j};                                      % propagate j to i Xform
        fh = Xup{j}' * fh;                                     % CRBA
        j = model.parent(j);                                   % next j
        J{i}(:,j) = Xj * S{j};                                 % jacobian (still in i's coordinates)
        H(i,j) = S{j}' * fh;                                   % CRBA
        H(j,i) = H(i,j);                                       % CRBA
    end
    

    J0{i} = Xr0{i} * J{i};                                     % jacobian (now it world coordinates)
end

                          
for i = 1:model.NLEGS      
    %disp(i);% loop through feet  
    j = model.b_foot(i);                                       % body containing foot
    [~,pf0{i}] = plux_2(model.Xfoot{i} * X0{j});                 % origin to foot translation, world coordinates
    Jf0{i} = [zeros(3,3) eye(3)] * Xr0{j} * model.Xfoot{i} * Xr0{j}' * J0{j}; % transform world jacobian to foot
    % alternatively:
    %Jf0{i} = [zeros(3,3) eye(3)] * Xr0{j} * model.Xfoot{i} * J{j}
    
    vf = Xr0{j} * model.Xfoot{i} * v{j};                      % foot velocity
    Jdfqd0{i} = [zeros(3,3) eye(3)] * Xr0{j} * model.Xfoot{i} * a_no_qdd{j} + skew(vf(1:3)) * vf(4:6);
    vvf{i} = vf;
end

