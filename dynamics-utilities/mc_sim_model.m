function [ model, jtype, Xtree, I, Xrot, Irot ] = mc_sim_model(  )
% This function creates the rigid body model struct for the mini cheetah
% for some reason code generation doesn't work properly with the graphics
% stuff, so there is a separate MiniCheetahGraphicsModel that should be
% used in the simulator.

% The code is based on Cheetah3FullRotorModel_mex, with some slight
% restructuring and potentially incorrect comments added...

% like the big cheetah, 0 is with the legs straight down
% x is forward, y is left, and z is up


%% Model Initialization
NLEGS = 4; % number of legs

model.gc.point = zeros(3,8 +2*NLEGS); % ground contact points
model.gc.body = zeros(1,8 + 2*NLEGS); % ground contact body numbers

NJJ = 3*NLEGS; % number of joints

Nb = 6; % body index (starts at 6, the floating base)

model.parent = zeros(1,6+NJJ); % parent contains the body number of the parent body
model.parent(1) = 0;           % parent of the first body is 0
model.gr     = zeros(1,6+NJJ); % gear ratios

% these paremeters exist for each body:
jtype  = repmat({'  '},6+NJJ,1);        % joint type (floating base, rotation...)
Xtree  = repmat({eye(6)},6+NJJ,1);      % 6x6 transform from parent to child
I      = repmat({zeros(6)},6+NJJ,1);    % spatial inertia of body, body frame
Xrot   = repmat({eye(6)},6+NJJ,1);      % parent to child's rotor
Irot   = repmat({zeros(6)},6+NJJ,1);    % spatial inertia of rotor, rotor frame


%% Mini Cheetah Geometry Constants: All Constants in this section
abadWidth = .049*2;   % distance along y axis between left/right abads
abadLength = 0.19*2;  % distance along x axis between front/back abads

abadRotorWidth = .049*2;
abadRotorLength = .125*2;

l0 = 0.062;  % abad link length
l0_r = 0.04; % abad to hip rotor length
l1 = 0.209;  % leg upper link length
l2 = 0.175;  % leg lower link length


body_mass = 3.3;      % mass of body, not including rotors/legs
body_CoM = [0 0 0]; % CoM of body, body coordinates

gear_ratio_knee = 9.33;
gear_ratio_hip = 6;
gear_ratio_abad = 6;


% abad locations
abad_location = repmat({zeros(3,1)},4,1);         % abad joint locations, body frame
abad_location{1} = [ abadLength -abadWidth 0]'/2; % front right: leg 1
abad_location{2} = [ abadLength  abadWidth 0]'/2; % front left:  leg 2
abad_location{3} = [-abadLength -abadWidth 0]'/2; % rear right:  leg 3
abad_location{4} = [-abadLength  abadWidth 0]'/2; % rear left:   leg 4

abad_r_location = repmat({zeros(3,1)},4,1);                   % abad rotor locations, body frame
abad_r_location{1} = [ abadRotorLength -abadRotorWidth 0]'/2; % front right: leg 1
abad_r_location{2} = [ abadRotorLength  abadRotorWidth 0]'/2; % front left:  leg 2
abad_r_location{3} = [-abadRotorLength -abadRotorWidth 0]'/2; % rear right:  leg 3
abad_r_location{4} = [-abadRotorLength  abadRotorWidth 0]'/2; % rear left:   leg 4

% --- inertia of base - using a box ---
body_inertia = 1e-6*[11253, 0, 0;     % Body inertia about COM
               0, 362030, 0;
               0, 0, 42673];

% --- collision bounding box for body ---
%%% Body bounding corners %%%
b_1 = [0.135, 0.05, 0.05];
b_2 = [0.135, 0.05, -0.05];
b_3 = [0.135, -0.05, 0.05];
b_4 = [0.135, -0.05, -0.05];
b_5 = [-0.135, 0.05, 0.05];
b_6 = [-0.135, 0.05, -0.05];
b_7 = [-0.135, -0.05, 0.05];
b_8 = [-0.135, -0.05, -0.05];

model.gc.point(:,1:8) = [b_1' b_2' b_3' b_4' b_5' b_6' b_7' b_8'];


% --- mass properties of rotor  ---
rotor_mass = 0.055; % not used.

rotor_inertia_z = 1e-6*[33, 0, 0;      % Rotor inertia, rotates about z-axis
                0, 33, 0;
                0, 0, 63];
% rotor which rotates around the x axis
rotor_inertia_x = ry(pi/2) * rotor_inertia_z * ry(pi/2)';
% rotor which rotates around the y axis
rotor_inertia_y = rx(pi/2) * rotor_inertia_z * rx(pi/2)';

% --- mass properties of abad link - using a box ---
abad_link_mass = 0.54;
abad_link_inertia = 1e-6*[381, 58, 0.45;        % Moment inertia about COM, aligned with ab/ad frame
              58, 560, 0.95;
              0.45, 0.95, 444];     
    

% the left and right abads are mirror images!
% the left ab/ad has the y axis point from center to hip
% the right ab/ad has the y axis point from hip to center
% so their CoM's relative to their origins are different
abad_CoM_left = [0 .036 0];
abad_CoM_right = [0 -.036 0];

% --- mass properies of hip (upper) link ---
hip_link_mass = 0.634;
hip_link_inertia = 1e-6*[1983, 245, 13;      % Moment of inertia about COM, aligned with hip frame
            245, 2103, 1.5;
            13, 1.5, 408];

% hip CoM is a mirror image again.
hip_CoM_left = [0 .016 -.02];
hip_CoM_right = [0 -.016 -.02];

% --- mass properties of knee (upper) link ---
knee_link_mass = 0.064;

klir = 1e-6*[6, 0, 0;        % Moment of inertia about COM, aligned with old knee frame
            0, 248, 0;
            0, 0, 245];
        
% need to rotate knee link inertia
knee_link_inertia = ry(pi/2) * klir * ry(pi/2)';
knee_CoM_left = [0 0 -.061];
knee_CoM_right = [0 0 -.061];



%% Geometry Stuff




%% Floating Base
model.parent(1:6) = 0;   % no parent
jtype{6} = 'Fb';         % floating base joint type
Xtree{6} = eye(6);       % identity transform for floating base
I{6}     = mcI(body_mass,body_CoM,body_inertia); % spatial inertia of body
Xrot{6}  = eye(6);       % identity for rotor transform - base has no rotors so this doesn't matter
Irot{6}  = zeros(6);     % rotor inertia is zero
model.gr(6) = 1;         % gear ratio - base has no rotors so this doesn't matter

% contact model for floating base
% all points belong to the base body...
model.gc.body(1:8) = [Nb, Nb, Nb, Nb, Nb, Nb, Nb, Nb];

gcInd = 9; % index for the ground contact points
% after the base, which has 8 points, it is 9

Nb = 6;          % body index (6 is floating base)
side_sign = -1;  % -1 for the first side (right)

%% LEGS
for i = 1:NLEGS
    %% Ab/Ad (hip roll)
    Nb = Nb + 1;                                 % increment body index
    model.parent(Nb) = 6;                        % ab/ad parent is 6, which is floating base
    model.gr(Nb) = gear_ratio_abad;                   % ab/ad gear ratio
    jtype{Nb}  = 'Rx';                           % rotate around x axis of floating base
    Xtree{Nb}  = plux(eye(3), abad_location{i}); % floating base to ab/ad
    Xrot{Nb}  =  plux(eye(3), abad_r_location{i}); % floating base to ab/ad rotor
    
    % pick left/right side ab/ad spatial inertia
    if(side_sign == -1)
        I{Nb} = mcI(abad_link_mass,abad_CoM_right,abad_link_inertia);
    else
        I{Nb} = mcI(abad_link_mass,abad_CoM_left,abad_link_inertia);
    end
    % rotor rotates around x axis for ab/ad
    Irot{Nb} = mcI(rotor_mass,[0 0 0],rotor_inertia_x);
    
    
    %% Hip (Hip Pitch)
    Nb = Nb + 1;
    model.parent(Nb) = Nb-1; % parent is ab/ad (hip gearbox rotated by ab/ad)
    model.gr(Nb) = gear_ratio_hip; 
    jtype{Nb}  = 'Ry'; % Hip rotates around y axis
    
    % Ab/ad to hip transform
    % rotate around z by pi so that positive hip angle moves the leg +x
    Xtree{Nb}  = plux(rz(pi), [0 side_sign*l0 0]');
    
    % Ab/ad to hip rotor transform is the same as ab/ad to hip transform
    Xrot{Nb} = plux(rz(pi), [0 side_sign*l0_r 0]');
    
    % pick left/right side hip spatial inertia
    if(side_sign == -1)
        I{Nb} = mcI(hip_link_mass,hip_CoM_right,hip_link_inertia);
    else
        I{Nb} = mcI(hip_link_mass,hip_CoM_left,hip_link_inertia);
    end
    
    % rotor rotates around y axis for hip
    Irot{Nb} = mcI(rotor_mass,[0 0 0],rotor_inertia_y);
    
    % ground contact knee
    model.gc.body(gcInd) = Nb;           % knee contact is on upper link
    model.gc.point(:,gcInd) = [0 0 -l1]; % hip origin to knee
    gcInd = gcInd + 1;                   % increment gcind
    
   %% Knee (Knee Pitch)
   Nb = Nb + 1;             % increment current joint number
   model.gr(Nb) = gear_ratio_knee; 
   model.parent(Nb) = Nb-1; % parent is hip
   jtype{Nb}  = 'Ry';       % knee rotates around y axis
   model.NB   = Nb;         % total number of joints can be set here

   Xtree{Nb}  = plux(eye(3), [0 0 -l1]'); % hip to knee transform
   Xrot{Nb}  = plux(eye(3), [0 0  0 ]');  % hip to knee rotor transform
   
    % pick left/right side hip spatial inertia
    if(side_sign == -1)
        I{Nb} = mcI(knee_link_mass,knee_CoM_right,knee_link_inertia);
    else
        I{Nb} = mcI(knee_link_mass,knee_CoM_left,knee_link_inertia);
    end
    model.gc.body(gcInd) = Nb;            % current body
    model.gc.point(:,gcInd) = [0 0 -l2]'; % add point at foot
    gcInd = gcInd + 1; % increment
    
    side_sign = -1 * side_sign; % flip side sign for the next leg.
end
end

function I = boxInertia(mass, x)
    I = (norm(x)^2*eye(3) - diag(x.^2))*mass/12;
end