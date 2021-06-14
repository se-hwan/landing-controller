function [  H,C,p0,J0,Jdqd,pf0,Jf0,Jdfqd0,vvf ] = dynamics_one_step( model, rotor_model, q, qd )
% compute dynamics... with ROTORS!

%% INITIALIZE QUANTITIES:
H = zeros(model.NB);
C = zeros(model.NB,1);
b_rot = repmat({zeros(6,1)},model.NB,1);
Ic = model.I;

a_grav = [0 0 0 0 0 -9.81]';          

%% LOOP 1 (through bodies, base to tip)
for i = 1:model.NB
    J{i} = zeros(6,model.NB);                           % initialize jacobian
    Jdqd{i} = zeros(6,1);                               % initialize jdqd
    [XJ, S{i}] = jcalc( model.jtype{i}, q(i) );         % joint xform and motion subspace
    vJ = S{i}*qd(i);                                    % joint spatial velocity
    Xup{i} = XJ * model.Xtree{i};                       % from up the tree xform
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


%% LOOP 2 (through rotors)
for k = 1:rotor_model.NR
    i = rotor_model.gamma(k);                                  % gearing constraint joint
    [XJ,~] = jcalc(model.jtype{i},q(i));                       % joint xform
    Xup_rotor = XJ * rotor_model.X_mu{k};                      % from up the tree xform
    H(i,i) = H(i,i) +... 
        rotor_model.gr(k)^2 * S{i}' * rotor_model.I{k} * S{i}; % diagonal mass matrix
    b_rot{i} = b_rot{i} +...
        rotor_model.gr(k) * rotor_model.I{k} * S{i};           % off-diagonal term
    vJ = rotor_model.gr(k) * S{i} * qd(i);                     % spatial velocity of rotor joint
    if(model.parent(i) == 0)
        vk = zeros(6,1);                     % spatial velocity of rotor
        ak = zeros(6,1);      % spatial acceleration of rotor (velocity product)
        fk = zeros(6,1); % spatial force on rotor (velocity product)
    
        C(i,1) = C(i,1) + rotor_model.gr(k) * S{i}' * fk;
    else 
        vk = Xup_rotor * v{model.parent(i)};                       % spatial velocity of rotor
        ak = Xup_rotor * avp{model.parent(i)} + crm(vk) * vJ;      % spatial acceleration of rotor (velocity product)
        fk = rotor_model.I{k} * ak + crf(vk) * rotor_model.I{k} * vJ; % spatial force on rotor (velocity product)
    
        C(i,1) = C(i,1) + rotor_model.gr(k) * S{i}' * fk;
        fvp{model.parent(i)} = fvp{model.parent(i)} + Xup_rotor' * fk;
        Ic{model.parent(i)} = Ic{model.parent(i)} + Xup_rotor' * rotor_model.I{k} * Xup_rotor;
    end
   
    
    if(model.parent(i) ~= rotor_model.mu(k))
        disp(['ERROR: parent(i) = ' num2str(model.parent(i)) ' but mu(k) = ' num2str(rotor_model.mu(k))]);
    end
end

%% LOOP 3 (through bodies, bottom to top)
for i = model.NB:-1:1  
    C(i,1) = C(i,1) + S{i}' * fvp{i};   
    if model.parent(i) ~= 0
        fvp{model.parent(i)} = fvp{model.parent(i)} + Xup{i}'*fvp{i};  
         Ic{model.parent(i)} = Ic{model.parent(i)} + Xup{i}'*Ic{i}*Xup{i}; 
    end
end

%% LOOP 4 (through bodies, CRBA)
for i = 1:model.NB
    Xj = eye(6);                                               % from j to i (right now j = i)
    J{i}(:,i) = S{i};                                          % jacobian with i to i
    Jdqd{i} = Xr0{i} * a_no_qdd{i};                            % Jdqd is spatial acceleration due to qd, rotated to world
    v0 = Xr0{i} * v{i};                                        % spatial velocity, rotated to world
    Jdqd{i}(4:6) = Jdqd{i}(4:6) + skew(v0(1:3)) * v0(4:6);     % spatial -> normal acceleration
    
    fh = Ic{i} * S{i};                                         % CRBA
    H(i,i) = S{i}' * fh;                                       % CRBA
    j = i;                  
    while model.parent(j) > 0                                  % loop through i's parents (up)
        Xj = Xj * Xup{j};                                      % propagate j to i Xform
        fh = Xup{j}' * fh;                                     % CRBA
        b_rot{i} = Xup{j}' * b_rot{i};                         % CRBA type xform for rotor off-diagonal rxn torque
        j = model.parent(j);                                   % next j
        J{i}(:,j) = Xj * S{j};                                 % jacobian (still in i's coordinates)
        H(i,j) = S{j}' * (fh + b_rot{i});                      % CRBA + rotor off-diagonal
        H(j,i) = H(i,j);                                       % CRBA
    end
    

    J0{i} = Xr0{i} * J{i};                                     % jacobian (now it world coordinates)
end

%% LOOP 5 (through feet, does kinematics)
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
end

