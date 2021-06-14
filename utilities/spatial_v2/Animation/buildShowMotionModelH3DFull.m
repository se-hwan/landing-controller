function model = buildShowMotionModelH3DFull(params)

%% Humanoid 3D
if(strcmp(params.model,'humanoid3D')) % 3D humanoid model
    % Initialize model struct:
    NLEGS = 2;
    NARMS = 2;
    N_GND_CONTACTS = 2*NLEGS;
    model.NB = 22 + 4*N_GND_CONTACTS;                % number of bodies
    model.NLEGS = NLEGS;                             % number of legs
    model.NARMS = NARMS;                             % number of arms
    model.N_GND_CONTACTS = N_GND_CONTACTS;           % number of ground contact points
    model.gc_parent_limb = zeros(1,N_GND_CONTACTS);  % which limb does gc_contact point correspond to
    model.gravity   = [0 0 -9.81];                   % gravity
    model.parent    = zeros(1,model.NB);             % parent body indices
    model.jtype     = repmat({'  '},model.NB,1);     % joint types
    model.Xtree     = repmat({eye(6)},model.NB,1);   % coordinate transforms
    model.I         = repmat({zeros(6)},model.NB,1); % spatial inertias
    model.Xfoot     = repmat({zeros(6)},NLEGS,1);    % feet
    model.b_foot    = zeros(1,NLEGS);
    model.Xtoe      = repmat({zeros(6)},NLEGS,1);    % toe
    model.b_toe     = zeros(1,NLEGS);
    model.Xheel     = repmat({zeros(6)},NLEGS,1);    % heel
    model.b_heel    = zeros(1,NLEGS);
    model.Xhand     = repmat({zeros(6)},NARMS,1);    % hand
    model.b_hand    = zeros(1,NARMS);
    q_arm = [0 0 -0.1];
    q_leg = [0 0 -0.4 0.77 -0.37];
    model.q_home = [zeros(1,6) repmat(q_leg,1,NLEGS) repmat(q_arm,1,NARMS)];
    
    nb = 0; % current body index
    
    % Base x translation (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;   % parent is previous
    model.jtype{nb}  = 'Px';     % prismatic x
    model.Xtree{nb}  = eye(6);   % on top of previous joint
    model.I{nb}      = zeros(6); % massless
    % Base y translation (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;   % parent is previous
    model.jtype{nb}  = 'Py';     % prismatic y
    model.Xtree{nb}  = eye(6);   % on top of previous joint
    model.I{nb}      = zeros(6); % massless
    % Base z translation (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Pz';
    model.Xtree{nb}  = eye(6);
    model.I{nb}      = zeros(6);
    % Base roll (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Rx';
    model.Xtree{nb}  = eye(6);
    model.I{nb}      = zeros(6);
    % Base pitch (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Ry';
    model.Xtree{nb}  = eye(6);
    model.I{nb}      = zeros(6);
    % Base yaw (body mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Rz';
    model.Xtree{nb}  = eye(6);
    model.I{nb}      = zeros(6);
    
    % Loop through legs
    side_sign = [1 1;1 -1;1 1];
    %leg_side = -1;
    nb_base = nb;
    for leg = 1:NLEGS
        % Hip Rz
        nb = nb + 1;
        model.parent(nb) = nb_base;
        model.jtype{nb}  = 'Rz';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.hipRzLocation);
        model.I{nb}      = zeros(6);
        
        % Hip Rx
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Rx';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.hipRxLocation);
        model.I{nb}      = zeros(6);
        
        % Hip Ry
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.hipRyLocation);
        model.I{nb}      = zeros(6);
        
        % Knee
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.kneeLocation);
        model.I{nb}      = zeros(6);
        
        % Ankle
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.ankleLocation);
        model.I{nb}      = zeros(6);
        
        % Foot (bonus)
        model.Xfoot{leg}   = plux(eye(3),[0 0 -params.footHeight+0.2]');
        model.b_foot(leg)  = nb;
        
        % Toe (bonus)
        model.Xtoe{leg}   = plux(eye(3),[params.footToeLength 0 -params.footHeight]');
        model.b_toe(leg)  = nb;
        model.gc_parent_limb(2*(NLEGS-1)+1) = leg;
        
        % Heel (bonus)
        model.Xheel{leg}   = plux(eye(3),[-params.footHeelLength 0 -params.footHeight]');
        model.b_heel(leg)  = nb;
        model.gc_parent_limb(2*(NLEGS-1)+2) = leg;
        
    end
    
    % Loop through arms
    %leg_side = -1;
    for arm = 1:NARMS
        % Shoulder Rx
        nb = nb + 1;
        model.parent(nb) = nb_base;
        model.jtype{nb}  = 'Rx';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.shoulderRxLocation);
        model.I{nb}      = zeros(6);
        
        % Shoulder Ry
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.shoulderRyLocation);
        model.I{nb}      = zeros(6);
        
        % Elbow
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.elbowLocation);
        model.I{nb}      = zeros(6);
        
        % Hand (bonus)
        model.Xhand{arm}   = plux(eye(3),[0 0 -params.lowerArmLength]');
        model.b_hand(arm)  = nb;
        
    end
    
    % Loop through the feet
    foot_force_parent = zeros(4,1);
    for leg = 1:model.NLEGS
        %%%% Toe
        % X direction
        nb = nb + 1;
        model.parent(nb) = 0;
        model.jtype{nb}  = 'Px';
        model.Xtree{nb}  = eye(6);
        model.I{nb}      = zeros(6,6);
        
        % Y direction
        nb = nb + 1;
        model.parent(nb) = nb-1;
        model.jtype{nb}  = 'Py';
        model.Xtree{nb}  = eye(6);
        model.I{nb}      = zeros(6,6);
        
        % Z direction
        nb = nb + 1;
        model.parent(nb) = nb-1;
        model.jtype{nb}  = 'Pz';
        model.Xtree{nb}  = eye(6);
        model.I{nb}      = zeros(6,6);
        foot_force_parent(2*leg-1) = nb;
        
        %%%% Heel
        % X direction
        nb = nb + 1;
        model.parent(nb) = 0;
        model.jtype{nb}  = 'Px';
        model.Xtree{nb}  = eye(6);
        model.I{nb}      = zeros(6,6);
        
        % Y direction
        nb = nb + 1;
        model.parent(nb) = nb-1;
        model.jtype{nb}  = 'Py';
        model.Xtree{nb}  = eye(6);
        model.I{nb}      = zeros(6,6);
        
        % Z direction
        nb = nb + 1;
        model.parent(nb) = nb-1;
        model.jtype{nb}  = 'Pz';
        model.Xtree{nb}  = eye(6);
        model.I{nb}      = zeros(6,6);
        foot_force_parent(2*leg) = nb;
    end
    
    % Loop Through Force "Arrows"
    for leg = 1:model.NLEGS
        % Toe
        nb = nb + 1;
        model.parent(nb) = foot_force_parent(2*leg-1);
        model.jtype{nb}  = 'Pz';
        model.Xtree{nb}  = eye(6);
        model.I{nb}      = zeros(6,6);
        % Heel
        nb = nb + 1;
        model.parent(nb) = foot_force_parent(2*leg);
        model.jtype{nb}  = 'Pz';
        model.Xtree{nb}  = eye(6);
        model.I{nb}      = zeros(6,6);
    end
    
    nb = 0;
    % Base x translation (no mass)
    nb = nb + 1;
    model.appearance.body{nb} = {};
    % Base y translation (no mass)
    nb = nb + 1;
    model.appearance.body{nb} = {};
    % Base z translation (no mass)
    nb = nb + 1;
    model.appearance.body{nb} = {};
    % Base roll (no mass)
    nb = nb + 1;
    model.appearance.body{nb} = {};
    % Base pitch (no mass)
    nb = nb + 1;
    model.appearance.body{nb} = {};
    % Base yaw (body mass)
    nb = nb + 1;
    model.appearance.body{nb} = ...
        {'colour',[1.0 1.0 1.0],...
        'box', [-0.5*params.torso_length -0.5*params.torso_width -0.0*params.torso_height;...
        0.5*params.torso_length 0.5*params.torso_width 1.0*params.torso_height]};
    
    side_sign = [1 1;1 -1;1 1];
    for leg = 1:model.NLEGS
        % Hip Rz
        nb = 7;
        cyl1 = [0 0 0.75*params.leg_rad;...
            0 0 -0.75*params.leg_rad];
        cyl2 = [0 0 0;...
            side_sign(:,leg)'.*params.hipRxLocation];
        model.appearance.body{nb+5*(leg-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 1.65*params.leg_rad,...
            'cyl', cyl2, params.leg_rad};
        
        % Hip Rx
        nb = 8;
        cyl1 = [0.5*params.leg_rad 0 0;...
            -0.5*params.leg_rad 0 0];
        cyl2 = [0 0 0;...
            side_sign(:,leg)'.*params.hipRyLocation];
        model.appearance.body{nb+5*(leg-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 1.45*params.leg_rad,...
            'cyl', cyl2, params.leg_rad};
        
        % Hip Ry
        nb = 9;
        cyl1 = [0 0.5*params.leg_rad 0;...
            0 -0.5*params.leg_rad 0];
        cyl2 = [0 0 0;...
            side_sign(:,leg)'.*params.kneeLocation];
        model.appearance.body{nb+5*(leg-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 1.25*params.leg_rad,...
            'cyl', cyl2, params.leg_rad};
        
        % Knee
        nb = 10;
        cyl1 = [0 0.5*params.leg_rad 0;...
            0 -0.5*params.leg_rad 0];
        cyl2 = [0 0 0;...
            side_sign(:,leg)'.*params.ankleLocation];
        model.appearance.body{nb+5*(leg-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 1.15*params.leg_rad,...
            'cyl', cyl2, params.leg_rad};
        
        % Ankle
        nb = 11;
        cyl1 = [0 0.5*params.leg_rad 0;...
            0 -0.5*params.leg_rad 0];
        cyl2 = [0 0 0;...
            side_sign(:,leg)'.*[params.footToeLength 0 -params.footHeight]];
        cyl3 = [0 0 0;...
            side_sign(:,leg)'.*[-params.footHeelLength 0 -params.footHeight]];
        model.appearance.body{nb+5*(leg-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 0.8*params.leg_rad,...
            'cyl', cyl2, 0.5*params.leg_rad,...
            'cyl', cyl3, 0.5*params.leg_rad};
    end
    
    for arm = 1:model.NARMS
        % Shoulder Rx
        nb = 17;
        cyl1 = [0.75*params.leg_rad 0 0;...
            -0.75*params.leg_rad 0 0];
        cyl2 = [0 0 0;...
            side_sign(:,arm)'.*params.shoulderRyLocation];
        model.appearance.body{nb+3*(arm-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 1.65*params.leg_rad,...
            'cyl', cyl2, params.leg_rad};
        
        % Shoulder Ry
        nb = 18;
        cyl1 = [0 0.75*params.leg_rad 0;...
            0 -0.75*params.leg_rad 0];
        cyl2 = [0 0 0;...
            side_sign(:,arm)'.*params.elbowLocation];
        model.appearance.body{nb+3*(arm-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 1.35*params.leg_rad,...
            'cyl', cyl2, params.leg_rad};
        
        % Elbow
        nb = 19;
        cyl1 = [0 0.75*params.leg_rad 0;...
            0 -0.75*params.leg_rad 0];
        cyl2 = [0 0 0;...
            side_sign(:,arm)'.*[0 0 -params.lowerArmLength]];
        model.appearance.body{nb+3*(arm-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 1.15*params.leg_rad,...
            'cyl', cyl2, params.leg_rad};
        
    end
    
    % Feet (R toe, R heel , L Toe, L Heel)
    model.appearance.body{25} = ...
        {'colour',[0.95 0.05 0.05],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{35} = ...
        {'colour',[0.95 0.05 0.05],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{28} = ...
        {'colour',[0 0.95 0.05],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{36} = ...
        {'colour',[0 0.95 0.05],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{31} = ...
        {'colour',[0 0.05 0.95],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{37} = ...
        {'colour',[0 0.05 0.95],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{34} = ...
        {'colour',[0.95 0.85 0.],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{38} = ...
        {'colour',[0.95 0.85 0.],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    
    % Ground
    x_start = -0.5;
    x_end = 2.0;
    y_width = [-1.0 1.0];
    samples = 100;
    x = linspace(x_start, x_end, samples);
    
    groundV = [];
    for k = 1:samples
        groundV = [groundV;x(k) y_width(1) getGroundHeight(params.obstacle, [x(k) y_width(1)])];
        groundV = [groundV;x(k) y_width(2) getGroundHeight(params.obstacle, [x(k) y_width(2)])];
    end
    
    N_triangles = size(groundV,1)-2;
    groundT = zeros(N_triangles, 3);
    for k = 1:N_triangles
        if (rem(k,2) == 0) % even
            groundT(k,:) = [k k+1 k+2];
        else
            groundT(k,:) = [k+1 k k+2];
        end
    end
    
    model.appearance.base = ...
        {'vertices', groundV, ...
        'triangles',groundT};
    
    
    disp(['Created robot model with ' num2str(model.NB) ' coordinates!']);
    
end