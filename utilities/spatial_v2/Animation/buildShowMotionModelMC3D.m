function model = buildShowMotionModelMC3D(params, model, ground_height)


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

% % Joints
% 7 - FR ab/ad
% 8 - FR hip
% 9 - FR knee
% 10 - FL ab/ad
% 11 - FL hip
% 12 - FL knee
% .
% .
% .

% Ab/Ad
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

% Ground
model.appearance.base = ...
    { 'box', [-1.5 -1.5 ground_height;...
    10.5 1.5 ground_height-0.01]};
