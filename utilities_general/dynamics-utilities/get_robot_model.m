function [ model ] = get_robot_model( params )
model.name = params.model;

% Get spatial_v2 style robot model for 2D quadruped or humanoid
%
% ************     2D Spatial Quadruped Model:           ************
% Coordinate system: x forward, y left, z up
% Bodies:
% 0   - fixed origin
% 1   - Px (translate horizontal) (massless)
% 2   - Pz (translate vertical)   (massless)
% 3   - Ry (pitch)                (base + ab/ad link, rotor)
% 4,5 - Front "leg"
% 6,7 - Rear "leg"
% Bonus - foot positions
%
% ************     3D Quadruped Model:           ************
% Coordinate system: x forward, y left, z up
% Bodies:
% 0   - fixed origin
% 1   - Px (translate fore/aft)   (massless)
% 2   - Py (translate lateral)    (massless)
% 3   - Pz (translate vertical)   (massless)
% 4   - Rx (roll)                 (massless)
% 5   - Ry (pitch)                (massless)
% 6   - Rz (yaw)                  (base + ab/ad link, rotor)
% 7,8,9 - Front right "leg"
% 10,11,12 - Front left "leg"
% 13,14,15 - Rear right leg
% 16,17,18 - Rear left "leg"
% Bonus - foot positions
%
% ************     Reduced Order Planar Quadruped Model:           ************
% Coordinate system: x forward, y up
% Bodies:
% 0   - fixed origin
% 1   - px (translate horizontal) (massless)
% 2   - pz (translate vertical)   (massless)
% 3   - r (pitch)                (massless)
% 4 - Front foot
% 5 - Rear Foot
%
% ************      2D Spatial Humanoid Model:           ************
% Coordinate system: x forward, y left, z up
% Bodies:
% 0   - fixed origin
% 1   - Px (translate horizontal) (massless)
% 2   - Pz (translate vertical)   (massless)
% 3   - Ry (pitch)                (base + ab/ad link, rotor)
% 4,5,6 - Legs
% 8,9 - Arms
% Bonus - foot and hand positions
%
% ************      3D Humanoid Model:           ************
% Coordinate system: x forward, y left, z up
% Bodies:
% 0   - fixed origin
% 1   - Px (translate fore/aft)   (massless)
% 2   - Py (translate lateral)    (massless)
% 3   - Pz (translate vertical)   (massless)
% 4   - Rx (roll)                 (massless)
% 5   - Ry (pitch)                (massless)
% 6   - Rz (yaw)                  (base + ab/ad link, rotor)
% 7,8,9,10,11 - Right Leg
% 12,13,14,15,16 - Left leg
% 17,18,19 - right arm
% 20,21,22 - left arm
% Bonus - foot (and hand?) positions

%% Quadruped
if(strcmp(params.model,'quad')) % 2D quadruped model in spatial coords
    % Various mass parameters and locations
    I_body = mcI(params.body_mass, [0 0 0], params.i_body_mult * boxInertia(params.body_mass,[params.body_length, params.body_width, params.body_height]));
    I_l1   = mcI(params.l1_mass, [0 0 -params.l1/2], boxInertia(params.l1_mass,[params.leg_rad*2 params.leg_rad*2 params.l1]));
    I_l2   = mcI(params.l2_mass, [0 0 -params.l2/2], boxInertia(params.l2_mass,[params.leg_rad*2 params.leg_rad*2 params.l2]));
    hip_x = [params.body_length -params.body_length]/2;
    NLEGS = 2;
    NARMS = 0;
    % Initialize model struct:
    model.NB = 7;                                  % number of bodies
    model.NLEGS = NLEGS;                           % number of legs
    model.NARMS = NARMS;                           % number of arms
    model.gravity = [0 0 -9.81];                   % gravity
    model.parent  = zeros(1,model.NB);             % parent body indices
    model.jtype   = repmat({'  '},model.NB,1);     % joint types
    model.Xtree   = repmat({eye(6)},model.NB,1);   % coordinate transforms
    model.I       = repmat({zeros(6)},model.NB,1); % spatial inertias
    model.Xfoot   = repmat({zeros(6)},NLEGS,1);    % feet
    model.b_foot  = zeros(1,NLEGS);
    nb = 0; % current body index
    
    % Base x translation (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;   % parent is previous
    model.jtype{nb}  = 'Px';     % prismatic x
    model.Xtree{nb}  = eye(6);   % on top of previous joint
    model.I{nb}      = zeros(6); % massless
    % Base z translation (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Pz';
    model.Xtree{nb}  = eye(6);
    model.I{nb}      = zeros(6);
    % Base pitch (body mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb} = 'Ry';
    model.Xtree{nb} = eye(6);
    model.I{nb}     = I_body;
    
    % Loop through legs (just two for now)
    nb_base = nb;
    for leg = 1:NLEGS
        % Hip
        nb = nb + 1;
        model.parent(nb) = nb_base;  % hip parent is base
        model.jtype{nb}  = 'Ry';     % rotate around y
        %                 flip so y+ is leg forward,   translate to front/back
        model.Xtree{nb}  = plux(rz(pi),[0 0 0]') * plux(eye(3),[hip_x(leg) 0 0]');
        model.I{nb}      = I_l1;
        % Knee
        nb = nb + 1;
        model.parent(nb) = nb - 1; % knee parent is hip
        model.jtype{nb}  = 'Ry'; % rotate around y
        model.Xtree{nb}  = plux(eye(3),[0 0 -params.l1]');
        model.I{nb}      = I_l2;
        % Foot (bonus)
        model.Xfoot{leg}  = plux(eye(3),[0 0 -params.l2]');
        model.b_foot(leg) = nb;
    end
    disp(['Created robot model with ' num2str(nb) ' coordinates!']);
    
    %% Quadruped 3D
elseif(strcmp(params.model,'quad3D')) % 3D quadruped model in spatial coords
    % Initialize model struct:
    model.NB = 18;                                    % number of bodies
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
    % Base yaw (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb} = 'Rz';
    model.Xtree{nb} = eye(6);
    model.I{nb}     = zeros(6);
    % Base pitch (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb} = 'Ry';
    model.Xtree{nb} = eye(6);
    model.I{nb}     = zeros(6);
    % Base roll (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb} = 'Rx';
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
        if leg_side > 0
            model.I{nb}      = params.abadInertia;
        else
            model.I{nb}      = flipAlongAxis(params.abadInertia,'Y');
        end
        
        % Hip
        nb = nb + 1;
        model.parent(nb) = nb-1;  % hip parent is base
        model.jtype{nb}  = 'Ry';     % rotate around y
        model.Xtree{nb}  = plux(rz(pi),[0 0 0]')*plux(eye(3),side_sign(:,leg)'.*params.hipLocation);
        if leg_side > 0
            model.I{nb}      = params.hipInertia;
        else
            model.I{nb}      = flipAlongAxis(params.hipInertia,'Y');
        end
        
        % Knee
        nb = nb + 1;
        model.parent(nb) = nb-1;  % hip parent is base
        model.jtype{nb}  = 'Ry';     % rotate around y
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.kneeLocation);
        if leg_side > 0
            model.I{nb}      = params.kneeInertia;
        else
            model.I{nb}      = flipAlongAxis(params.kneeInertia,'Y');
        end
        
        % Foot (bonus)
        model.Xfoot{leg}  = plux(eye(3),side_sign(:,leg)'.*params.footLocation);
        model.b_foot(leg) = nb;
        model.gc_parent_limb(leg) = leg;
        
        leg_side = -1 * leg_side;
    end
    
    % Actuators
    model.gr = [params.abadGearRatio;params.hipGearRatio;params.kneeGearRatio];
    model.kt = [params.motorKT;params.motorKT;params.motorKT];
    model.Rm = [params.motorR;params.motorR;params.motorR];
    tauMax = model.gr.*[params.motorTauMax;params.motorTauMax;params.motorTauMax];
    model.tauMax = repmat(tauMax,4,1);
    model.batteryV = params.batteryV;
    
    disp(['Created robot model with ' num2str(nb) ' coordinates!']);
    
    %% Humanoid
elseif(strcmp(params.model,'humanoid')) % 2D humanoid model in full spatial coordinates
    % Various mass parameters and locations
    I_body = mcI(params.body_mass, [0 0 0], boxInertia(params.body_mass,[params.body_length, params.body_width, params.body_height]));
    
    I_shoulder  = mcI(params.shoulderMassProperties(1),...
        [0 0 -params.l_upper_arm/2],...
        boxInertia(params.shoulderMassProperties(1),[params.upper_arm_rad*2 params.upper_arm_rad*2 params.l_upper_arm]));
    I_elbow   = mcI(params.elbowMassProperties(1),...
        [0 0 -params.l_forearm/2],...
        boxInertia(params.elbowMassProperties(1),[params.forearm_rad*2 params.forearm_rad*2 params.l_forearm]));
    
    I_hip   = I_shoulder;
    I_knee  = I_elbow;
    I_ankle = I_elbow;
    
    % Initialize model struct:
    NLEGS = 1;
    NARMS = 1;
    model.NB = 8;                                  % number of bodies
    model.NLEGS = NLEGS;                           % number of legs
    model.NARMS = NARMS;                           % number of arms
    model.gravity   = [0 0 -9.81];                   % gravity
    model.parent    = zeros(1,model.NB);             % parent body indices
    model.jtype     = repmat({'  '},model.NB,1);     % joint types
    model.Xtree     = repmat({eye(6)},model.NB,1);   % coordinate transforms
    model.I         = repmat({zeros(6)},model.NB,1); % spatial inertias
    model.Xfoot     = repmat({zeros(6)},NLEGS,1);    % feet
    model.b_foot    = zeros(1,NLEGS);
    model.Xhand     = repmat({zeros(6)},NARMS,1);    % hand
    model.b_hand    = zeros(1,NARMS);
    nb = 0; % current body index
    
    % Fixed robot based for animation
    %     model.appearance.base = ...
    %         { 'box', [-0.05*params.body_length -0.05*params.body_width -0.05*params.body_height;...
    %         0.05*params.body_length 0.05*params.body_width 0.05*params.body_height]};
    
    % Base x translation (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;   % parent is previous
    model.jtype{nb}  = 'Px';     % prismatic x
    model.Xtree{nb}  = eye(6);   % on top of previous joint
    model.I{nb}      = zeros(6); % massless
    % Base z translation (no mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Pz';
    model.Xtree{nb}  = eye(6);
    model.I{nb}      = zeros(6);
    % Base pitch (body mass)
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb} = 'Ry';
    model.Xtree{nb} = eye(6);
    model.I{nb}     = I_body;
    
    % Loop through arms
    nb_base = nb;
    for arm = 1:NARMS
        % Shoulder
        nb = nb + 1;
        model.parent(nb) = nb_base;  % hip parent is base
        model.jtype{nb}  = 'Ry';     % rotate around y
        model.Xtree{nb}  = plux(eye(3),params.shoulder_location');
        model.I{nb}      = I_shoulder;
        % Elbow
        nb = nb + 1;
        model.parent(nb) = nb - 1; % knee parent is hip
        model.jtype{nb}  = 'Ry'; % rotate around y
        model.Xtree{nb}  = plux(eye(3),params.elbow_location');
        model.I{nb}      = I_elbow;
        % Hand (bonus)
        model.Xhand{arm}   = plux(eye(3),[0 0 -params.l_forearm]');
        model.b_hand(arm)  = nb;
    end
    
    % Loop through legs (only one in 2D model)
    for leg = 1:NLEGS
        % Hip
        nb = nb + 1;
        model.parent(nb) = nb_base;  % hip parent is base
        model.jtype{nb}  = 'Ry';     % rotate around y
        model.Xtree{nb}  = plux(eye(3),params.hip_location');
        model.I{nb}      = I_hip;
        % Knee
        nb = nb + 1;
        model.parent(nb) = nb - 1; % knee parent is hip
        model.jtype{nb}  = 'Ry'; % rotate around y
        model.Xtree{nb}  = plux(eye(3),params.knee_location');
        model.I{nb}      = I_knee;
        % Ankle
        nb = nb + 1;
        model.parent(nb) = nb - 1; % ankle parent is knee
        model.jtype{nb}  = 'Ry'; % rotate around y
        model.Xtree{nb}  = plux(eye(3),params.ankle_location');
        model.I{nb}      = I_ankle;
        
        % Foot (bonus)
        model.Xfoot{leg}   = plux(eye(3),[0 0 -params.foot_height]'); % contact point below ankle
        model.b_foot(leg)  = nb;
    end
    
    disp(['Created robot model with ' num2str(nb) ' coordinates!']);
    
    %% Humanoid 3D
elseif(strcmp(params.model,'humanoid3D')) % 3D humanoid model
    % Various mass parameters and locations
    I_body = mcI(params.torso_mass, params.torso_COM, boxInertia(params.torso_mass,[params.torso_length, params.torso_width, params.torso_height]));
    I_shoulderRx = mcI(params.shoulderRxMass,params.shoulderRxCOM,reshape(params.shoulderRxRotationalInertia,3,3));
    I_shoulderRy = mcI(params.shoulderRyMass,params.shoulderRyCOM,reshape(params.shoulderRyRotationalInertia,3,3));
    I_elbow = mcI(params.elbowMass,params.elbowCOM,reshape(params.elbowRotationalInertia,3,3));
    I_hipRz = mcI(params.hipRzMass,params.hipRzCOM,reshape(params.hipRzRotationalInertia,3,3));
    I_hipRx = mcI(params.hipRxMass,params.hipRxCOM,reshape(params.hipRxRotationalInertia,3,3));
    I_hipRy = mcI(params.hipRyMass,params.hipRyCOM,reshape(params.hipRyRotationalInertia,3,3));
    I_knee = mcI(params.kneeMass,params.kneeCOM,reshape(params.kneeRotationalInertia,3,3));
    I_ankle = mcI(params.ankleMass,params.ankleCOM,reshape(params.ankleRotationalInertia,3,3));
    
    %     I_shoulderRx  = massPropToSpatialInertia(params.shoulderRxMassProperties);
    %     I_shoulderRy  = massPropToSpatialInertia(params.shoulderRyMassProperties);
    %     I_elbow  = massPropToSpatialInertia(params.elbowMassProperties);
    %
    %     I_hipRz   = I_shoulderRx;
    %     I_hipRx   = I_shoulderRy;
    %     I_hipRy   = I_shoulderRy;
    %     I_knee  = I_elbow;
    %     I_ankle = I_elbow;
    
    % Initialize model struct:
    NLEGS = 2;
    NARMS = 2;
    N_GND_CONTACTS = 2*NLEGS;
    model.NB = 22;                                   % number of bodies
    model.NLEGS = NLEGS;                             % number of legs
    model.NARMS = NARMS;                             % number of arms
    model.N_GND_CONTACTS = N_GND_CONTACTS;           % number of ground contact points
    model.gc_parent_limb = zeros(1,N_GND_CONTACTS);          % which limb does gc_contact point correspond to
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
    model.I{nb}      = I_body;
    
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
        if leg == 2
            model.I{nb}      = I_hipRz;
        else
            model.I{nb}      = flipAlongAxis(I_hipRz,'Y');
        end
        
        % Hip Rx
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Rx';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.hipRxLocation);
        if leg == 2
            model.I{nb}      = I_hipRx;
        else
            model.I{nb}      = flipAlongAxis(I_hipRx,'Y');
        end
        
        % Hip Ry
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.hipRyLocation);
        if leg == 2
            model.I{nb}      = I_hipRy;
        else
            model.I{nb}      = flipAlongAxis(I_hipRy,'Y');
        end
        
        % Knee
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.kneeLocation);
        if leg == 2
            model.I{nb}      = I_knee;
        else
            model.I{nb}      = flipAlongAxis(I_knee,'Y');
        end
        
        % Ankle
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.ankleLocation);
        if leg == 2
            model.I{nb}      = I_ankle;
        else
            model.I{nb}      = flipAlongAxis(I_ankle,'Y');
        end
        
        % Foot (bonus)
        model.Xfoot{leg}   = plux(eye(3),[0 0 -params.footHeight+0.2]');
        model.b_foot(leg)  = nb;
        
        % Toe (bonus)
        model.Xtoe{leg}   = plux(eye(3),[params.footToeLength 0 -params.footHeight]');
        model.b_toe(leg)  = nb;
        model.gc_parent_limb(2*leg-1) = leg;
        
        % Heel (bonus)
        model.Xheel{leg}   = plux(eye(3),[-params.footHeelLength 0 -params.footHeight]');
        model.b_heel(leg)  = nb;
        model.gc_parent_limb(2*leg) = leg;
        
        %leg_side = -1 * leg_side;
    end
    
    % Loop through arms
    %leg_side = -1;
    for arm = 1:NARMS
        % Shoulder Rx
        nb = nb + 1;
        model.parent(nb) = nb_base;
        model.jtype{nb}  = 'Rx';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.shoulderRxLocation);
        if arm == 2
            model.I{nb}      = I_shoulderRx;
        else
            model.I{nb}      = flipAlongAxis(I_shoulderRx,'Y');
        end
        
        % Shoulder Ry
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.shoulderRyLocation);
        if arm == 2
            model.I{nb}      = I_shoulderRy;
        else
            model.I{nb}      = flipAlongAxis(I_shoulderRy,'Y');
        end
        
        % Elbow
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.elbowLocation);
        if arm == 2
            model.I{nb}      = I_elbow;
        else
            model.I{nb}      = flipAlongAxis(I_elbow,'Y');
        end
        
        % Hand (bonus)
        model.Xhand{arm}   = plux(eye(3),[0 0 -params.lowerArmLength]');
        model.b_hand(arm)  = nb;
        
        %leg_side = -1 * leg_side;
    end
    
    % Actuators (leg + arm)
    model.gr = [params.hipGearRatio;params.hipGearRatio;...
        params.hipGearRatio;params.kneeGearRatio;params.ankleGearRatio;...
        params.shoulderGearRatio;params.shoulderGearRatio;params.elbowGearRatio];
    model.kt = [params.smallMotorKT;params.smallMotorKT;...
        params.largeMotorKT;params.largeMotorKT;params.smallMotorKT;...
        params.smallMotorKT;params.smallMotorKT;params.smallMotorKT];
    model.Rm = [params.smallMotorR;params.smallMotorR;...
        params.largeMotorR;params.largeMotorR;params.smallMotorR;...
        params.smallMotorR;params.smallMotorR;params.smallMotorR];
    tauMax = model.gr.*[params.smallMotorTauMax;params.smallMotorTauMax;...
        params.largeMotorTauMax;params.largeMotorTauMax;params.smallMotorTauMax;...
        params.smallMotorTauMax;params.smallMotorTauMax;params.smallMotorTauMax];
    model.tauMax = [tauMax;tauMax];
    model.batteryV = params.batteryV;
    
    % index selection
    %floating base
    model.indx.x = [1 4];
    model.indx.y = [2 5];
    model.indx.z = [3 6];
    model.indx.rpy = [4 5 6];
    
    %legs
    model.indx.hipz = [7 12];
    model.indx.hipx = [8 13];
    model.indx.hipy = [9 14];
    model.indx.knee = [10 15];
    model.indx.ankle = [11 16];
    %arms
    model.indx.shoulderx = [17 20];
    model.indx.shouldery = [18 21];
    model.indx.elbow = [19 22];
    
    
    model.indx.leg = [model.indx.hipx model.indx.hipy model.indx.hipz model.indx.knee model.indx.ankle];
    model.indx.arm = [model.indx.shoulderx model.indx.shouldery model.indx.elbow model.indx.rpy];
    
    
    
    disp(['Created robot model with ' num2str(nb) ' coordinates!']);
    
    %% Humanoid* 3D
elseif(strcmp(params.model,'humanoid*3D')) % 3D humanoid model
    % Various mass parameters and locations
    I_body = mcI(params.torso_mass, params.torso_COM, boxInertia(params.torso_mass,[params.torso_length, params.torso_width, params.torso_height]));
    I_shoulderRz = mcI(params.shoulderRzMass,params.shoulderRzCOM,reshape(params.shoulderRzRotationalInertia,3,3));
    I_shoulderRx = mcI(params.shoulderRxMass,params.shoulderRxCOM,reshape(params.shoulderRxRotationalInertia,3,3));
    I_shoulderRy = mcI(params.shoulderRyMass,params.shoulderRyCOM,reshape(params.shoulderRyRotationalInertia,3,3));
    I_elbow = mcI(params.elbowMass,params.elbowCOM,reshape(params.elbowRotationalInertia,3,3));
    I_hipRz = mcI(params.hipRzMass,params.hipRzCOM,reshape(params.hipRzRotationalInertia,3,3));
    I_hipRx = mcI(params.hipRxMass,params.hipRxCOM,reshape(params.hipRxRotationalInertia,3,3));
    I_hipRy = mcI(params.hipRyMass,params.hipRyCOM,reshape(params.hipRyRotationalInertia,3,3));
    I_knee = mcI(params.kneeMass,params.kneeCOM,reshape(params.kneeRotationalInertia,3,3));
    I_ankle = mcI(params.ankleMass,params.ankleCOM,reshape(params.ankleRotationalInertia,3,3));
    
    % Initialize model struct:
    NLEGS = 2;
    NARMS = 2;
    model.NB = 24;                                  % number of bodies
    model.NLEGS = NLEGS;                           % number of legs
    model.NARMS = NARMS;                           % number of arms
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
    model.I{nb}      = I_body;
    
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
        if leg == 2
            model.I{nb}      = I_hipRz;
        else
            model.I{nb}      = flipAlongAxis(I_hipRz,'Y');
        end
        
        % Hip Rx
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Rx';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.hipRxLocation);
        if leg == 2
            model.I{nb}      = I_hipRx;
        else
            model.I{nb}      = flipAlongAxis(I_hipRx,'Y');
        end
        
        % Hip Ry
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.hipRyLocation);
        if leg == 2
            model.I{nb}      = I_hipRy;
        else
            model.I{nb}      = flipAlongAxis(I_hipRy,'Y');
        end
        
        % Knee
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.kneeLocation);
        if leg == 2
            model.I{nb}      = I_knee;
        else
            model.I{nb}      = flipAlongAxis(I_knee,'Y');
        end
        
        % Ankle
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.ankleLocation);
        if leg == 2
            model.I{nb}      = I_ankle;
        else
            model.I{nb}      = flipAlongAxis(I_ankle,'Y');
        end
        
        % Foot (bonus)
        model.Xfoot{leg}   = plux(eye(3),[0 0 -params.footHeight+0.2]');
        model.b_foot(leg)  = nb;
        
        % Toe (bonus)
        model.Xtoe{leg}   = plux(eye(3),[params.footToeLength 0 -params.footHeight]');
        model.b_toe(leg)  = nb;
        
        % Heel (bonus)
        model.Xheel{leg}   = plux(eye(3),[-params.footHeelLength 0 -params.footHeight]');
        model.b_heel(leg)  = nb;
        
        %leg_side = -1 * leg_side;
    end
    
    % Loop through arms
    %leg_side = -1;
    for arm = 1:NARMS
        % Shoulder Rz
        nb = nb + 1;
        model.parent(nb) = nb_base;
        model.jtype{nb}  = 'Rz';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.shoulderRzLocation);
        if arm == 2
            model.I{nb}      = I_shoulderRz;
        else
            model.I{nb}      = flipAlongAxis(I_shoulderRz,'Y');
        end
        
        % Shoulder Rx
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Rx';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.shoulderRxLocation);
        if arm == 2
            model.I{nb}      = I_shoulderRx;
        else
            model.I{nb}      = flipAlongAxis(I_shoulderRx,'Y');
        end
        
        % Shoulder Ry
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.shoulderRyLocation);
        if arm == 2
            model.I{nb}      = I_shoulderRy;
        else
            model.I{nb}      = flipAlongAxis(I_shoulderRy,'Y');
        end
        
        % Elbow
        nb = nb + 1;
        model.parent(nb) = nb - 1;
        model.jtype{nb}  = 'Ry';
        model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.elbowLocation);
        if arm == 2
            model.I{nb}      = I_elbow;
        else
            model.I{nb}      = flipAlongAxis(I_elbow,'Y');
        end
        
        % Hand (bonus)
        model.Xhand{arm}   = plux(eye(3),[0 0 -params.lowerArmLength]');
        model.b_hand(arm)  = nb;
        
        %leg_side = -1 * leg_side;
    end
    
    % Actuators (leg + arm)
    model.gr = [params.hipGearRatio;params.hipGearRatio;...
        params.hipGearRatio;params.kneeGearRatio;params.ankleGearRatio;...
        params.shoulderGearRatio;params.shoulderGearRatio;params.shoulderGearRatio;params.elbowGearRatio];
    model.kt = [params.smallMotorKT;params.smallMotorKT;...
        params.largeMotorKT;params.largeMotorKT;params.smallMotorKT;...
        params.smallMotorKT;params.smallMotorKT;params.smallMotorKT;params.smallMotorKT];
    model.Rm = [params.smallMotorR;params.smallMotorR;...
        params.largeMotorR;params.largeMotorR;params.smallMotorR;...
        params.smallMotorR;params.smallMotorR;params.smallMotorR;params.smallMotorR];
    tauMax = model.gr.*[params.smallMotorTauMax;params.smallMotorTauMax;...
        params.largeMotorTauMax;params.largeMotorTauMax;params.smallMotorTauMax;...
        params.smallMotorTauMax;params.smallMotorTauMax;params.smallMotorTauMax;params.smallMotorTauMax];
    model.tauMax = [tauMax;tauMax];
    model.batteryV = params.batteryV;
    
    disp(['Created robot model with ' num2str(nb) ' coordinates!']);
    
    
else
    error(['[get_robot_model] Invalid model type specified. Use ``quad"',...
        'for 2D quadruped model or ``humanoid" for 2D humanoid model'])
    
end


end % end of get_robot_model


% from Pat's model of cheetah 3
function I = boxInertia(mass, x)
I = (norm(x)^2*eye(3) - diag(x.^2))*mass/12;
end

function I = massPropToSpatialInertia(M)
I(1,1) = M(5);
I(1,2) = M(10);
I(1,3) = M(9);
I(2,1) = M(10);
I(2,2) = M(6);
I(2,3) = M(8);
I(3,1) = M(9);
I(3,2) = M(8);
I(3,3) = M(7);
cSkew = skew([M(2),M(3),M(4)]);
I(1:3,4:6) = cSkew;
I(4:6,1:3) = cSkew';
I(4:6,4:6) = M(1) * eye(3);
end

function I = flipAlongAxis(I_in, axis)
h = skew_spatial(I_in(1:3,4:6));
Ibar = I_in(1:3,1:3);
m = I_in(6,6);

if strcmp(class(I_in),'casadi.MX') % make sure this function works for optimization parameters
    P = casadi.MX.zeros(4,4);
    I = casadi.MX.eye(6);
else
    P = zeros(4,4);
    I = eye(6);
end

P(1:3,1:3) = 0.5*trace(Ibar)*eye(3) - Ibar;
P(1:3,4) = h;
P(4,1:3) = h';
P(4,4) = m;

X = eye(4);
if (axis == 'X')
    X(1, 1) = -1;
elseif (axis == 'Y')
    X(2, 2) = -1;
elseif (axis == 'Z')
    X(3, 3) = -1;
end
P = X * P * X;

% I = eye(6);
m = P(4,4);
h = P(1:3,4);
E = P(1:3,1:3);
Ibar = trace(E) * eye(3) - E;
I(1:3,1:3) = Ibar;
I(1:3,4:6) = skew(h);
I(4:6,1:3) = skew(h)';
I(4:6,4:6) = m * eye(3);
end
