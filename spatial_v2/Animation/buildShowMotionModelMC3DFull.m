function model = buildShowMotionModelMC3DFull(params)

%% Quadruped 3D
if(strcmp(params.model,'quad3D')) % 3D quadruped model in spatial coords
    % Initialize model struct:
    model.NB = 18 + 3*4 + 4;                          % number of bodies
    model.NLEGS = 4;                                  % number of legs
    model.N_GND_CONTACTS = model.NLEGS;               % number of ground contacts
    model.NARMS = 0;                                  % number of arms
    model.gravity = [0 0 -9.81];                      % gravity
    model.parent  = zeros(1,model.NB);                % parent body indices
    model.jtype   = repmat({'  '},model.NB,1);        % joint types
    model.Xtree   = repmat({eye(6)},model.NB,1);      % coordinate transforms
    model.I       = repmat({zeros(6)},model.NB,1);    % spatial inertias
    model.Xfoot   = repmat({zeros(6)},model.NLEGS,1); % feet
    model.b_foot  = zeros(1,model.NLEGS);
    model.gc_parent_limb = zeros(1,model.N_GND_CONTACTS); % which limb does gc_contact point correspond to
    q_leg = [0 -1.45 2.65];
    model.q_home  = [zeros(1,6) repmat(q_leg,1,model.NLEGS)];
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
    model.jtype{nb}  = 'Py';     % prismatic x
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
    model.jtype{nb} = 'Rx';
    model.Xtree{nb} = eye(6);
    model.I{nb}     = zeros(6);
    % Base pitch (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb} = 'Ry';
    model.Xtree{nb} = eye(6);
    model.I{nb}     = zeros(6);
    % Base yaw (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb} = 'Rz';
    model.Xtree{nb} = eye(6);
    model.I{nb}     = params.bodyInertia;
    
    % Loop through legs
    nb_base = nb;
    side_sign = [1 1 -1 -1;-1 1 -1 1;1 1 1 1];
    leg_side = -1;
    for leg = 1:model.NLEGS
        % Ab/Ad
        nb = nb + 1;
        model.parent(nb) = nb_base;  % hip parent is base
        model.jtype{nb}  = 'Rx';     % rotate around y
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.abadLocation);
        
        % Hip
        nb = nb + 1;
        model.parent(nb) = nb-1;  % hip parent is base
        model.jtype{nb}  = 'Ry';     % rotate around y
        model.Xtree{nb}  = plux(rz(pi),[0 0 0]')*plux(eye(3),side_sign(:,leg)'.*params.hipLocation);
        
        % Knee
        nb = nb + 1;
        model.parent(nb) = nb-1;  % hip parent is base
        model.jtype{nb}  = 'Ry';     % rotate around y
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.kneeLocation);
        
        % Foot (bonus)
        model.Xfoot{leg}  = plux(eye(3),side_sign(:,leg)'.*params.footLocation);
        model.b_foot(leg) = nb;
        model.gc_parent_limb(leg) = leg;
        
        leg_side = -1 * leg_side;
    end
    
    % Loop through the feet
    foot_force_parent = zeros(4,1);
    for leg = 1:model.NLEGS
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
        foot_force_parent(leg) = nb;
    end
    
    % Loop Through Force "Arrows"
    for leg = 1:model.NLEGS
        nb = nb + 1;
        model.parent(nb) = foot_force_parent(leg);
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
        'box', [-0.5*params.bodyLength -0.5*params.bodyWidth -0.5*params.bodyHeight;...
        0.5*params.bodyLength 0.5*params.bodyWidth 0.5*params.bodyHeight]};
        
    side_sign = [1 1 -1 -1;-1 1 -1 1;1 1 1 1];
    nb = 7;
    for leg = 1:4
        % Ab/ad
        nb = 7;
        cyl1 = [0.75*params.leg_rad 0 0;...
            -0.75*params.leg_rad 0 0];
        cyl2 = [0 0 0;...
            side_sign(:,leg)'.*params.hipLocation];
        model.appearance.body{nb+3*(leg-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 1.65*params.leg_rad,...
            'cyl', cyl2, params.leg_rad};
        
        % Hip
        nb = 8;
        cyl1 = [0 0.5*params.leg_rad 0;...
            0 -0.5*params.leg_rad 0];
        cyl2 = [0 0 0;...
            side_sign(:,leg)'.*params.kneeLocation];
        model.appearance.body{nb+3*(leg-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 1.45*params.leg_rad,...
            'cyl', cyl2, params.leg_rad};
        
        % Knee
        nb = 9;
        cyl1 = [0 0.5*params.leg_rad 0;...
            0 -0.5*params.leg_rad 0];
        cyl2 = [0 0 0;...
            side_sign(:,leg)'.*params.footLocation];
        model.appearance.body{nb+3*(leg-1)} = ...
            {'colour',[0.8 0.2 0.3],...
            'cyl', cyl1, 1.25*params.leg_rad,...
            'cyl', cyl2, params.leg_rad};
    end
    
    % Feet (FR, FL , BR, BL)
    model.appearance.body{21} = ...
        {'colour',[0.95 0.05 0.05],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{31} = ...
        {'colour',[0.95 0.05 0.05],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{24} = ...
        {'colour',[0 0.95 0.05],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{32} = ...
        {'colour',[0 0.95 0.05],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{27} = ...
        {'colour',[0 0.05 0.95],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{33} = ...
        {'colour',[0 0.05 0.95],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{30} = ...
        {'colour',[0.95 0.85 0.],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    model.appearance.body{34} = ...
        {'colour',[0.95 0.85 0.],...
        'sphere', [0 0 0], 1.65*params.leg_rad};
    
    % Ground
    try
        x_start = -1.0;
        x_end = 5.0;
        y_width = [-2.0 2.0];
        samples = 300;
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
    catch
        model.appearance.base = {'tiles', [-10 20;-6 6;0 0], 0.4};
    end
    
    % Camera
    model.camera.body = 6;
    model.camera.direction = [0.05 -0.6 0.2];
    model.camera.zoom = 0.65;
    
    disp(['Created animation model with ' num2str(model.NB) ' coordinates!']);
    
else
    error('Not ready for model other than quadruped for SRMB show motion')
    
end