function [ model ] = get_humanoid_arm_model()
% Get spatial_v2 style robot model for leg of humanoid robot
% Care only about kinematics
% All symbolic
% Right arn only

% 3D humanoid model
% Various mass parameters and locations
Id = eye(6);

% Initialize model struct:
model.NB = 3;                                  % number of bodies
model.gravity   = [0 0 -9.81];                   % gravity
model.parent    = zeros(1,model.NB);             % parent body indices
model.jtype     = repmat({'  '},model.NB,1);     % joint types
model.Xtree     = repmat({eye(6)},model.NB,1);   % coordinate transforms
model.I         = repmat({zeros(6)},model.NB,1); % spatial inertias
model.Xhand     = zeros(6);                      % feet
model.b_hand    = 0;
nb = 0; % current body index

hip_y = sym('hip_y','real');
knee_z = sym('knee_z','real');
foot_z = sym('foot_z','real');

% Loop through legs
nb_base = nb;
% Hip Rx
nb = nb + 1;
model.parent(nb) = nb_base;
model.jtype{nb}  = 'Rx';
model.Xtree{nb}  = eye(6);
model.I{nb}      = Id;

% Hip Ry
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Ry';
model.Xtree{nb}  = plux(eye(3),[0;hip_y;0]);
model.I{nb}      = Id;

% Knee
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Ry';
model.Xtree{nb}  = plux(eye(3),[0;0;knee_z]);
model.I{nb}      = Id;

% Foot (bonus)
model.Xfoot   = plux(eye(3),[0;0;foot_z]);
model.b_foot  = nb;

