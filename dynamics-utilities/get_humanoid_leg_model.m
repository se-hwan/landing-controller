function [ model ] = get_humanoid_leg_model(hipRx,hipRy,knee,ankle,foot)
% Get spatial_v2 style robot model for leg of humanoid robot
% Care only about kinematics
% All symbolic
% Right leg only

% 3D humanoid model
% Various mass parameters and locations
Id = eye(6);

% Initialize model struct:
model.NB = 5;                                  % number of bodies
model.gravity   = [0 0 -9.81];                   % gravity
model.parent    = zeros(1,model.NB);             % parent body indices
model.jtype     = repmat({'  '},model.NB,1);     % joint types
model.Xtree     = repmat({eye(6)},model.NB,1);   % coordinate transforms
model.I         = repmat({zeros(6)},model.NB,1); % spatial inertias
model.Xfoot     = zeros(6);                      % feet
model.b_foot    = 0;
nb = 0; % current body index

% Loop through legs
nb_base = nb;
% Hip Rz
nb = nb + 1;
model.parent(nb) = nb_base;
model.jtype{nb}  = 'Rz';
model.Xtree{nb}  = eye(6);
model.I{nb}      = Id;

% Hip Rx
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Rx';
model.Xtree{nb}  = plux(eye(3),hipRx);
model.I{nb}      = Id;

% Hip Ry
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Ry';
model.Xtree{nb}  = plux(eye(3),hipRy);
model.I{nb}      = Id;

% Knee
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Ry';
model.Xtree{nb}  = plux(eye(3),knee);
model.I{nb}      = Id;

% Ankle
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Ry';
model.Xtree{nb}  = plux(eye(3),ankle);
model.I{nb}      = Id;

% Foot (bonus)
model.Xfoot   = plux(eye(3),foot);
model.b_foot  = nb;

