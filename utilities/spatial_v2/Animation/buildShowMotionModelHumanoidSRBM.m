function model = buildShowMotionModelHumanoidSRBM(params, model, ground_height)


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
    'box', [-1.5*params.torso_length -0.9*params.torso_width -0.75*params.torso_height;...
    1.5*params.torso_length 0.9*params.torso_width 0.75*params.torso_height]};

side_sign = [1 1;1 -1;1 1];
for leg = 1:model.NLEGS
    % Hip Rz
    nb = 7;
    cyl1 = [0 0 0.75*params.leg_rad;...
        0 0 -0.75*params.leg_rad];
    cyl2 = [0 0 0;...
        side_sign(:,leg)'.*params.hipRxLocation];
    model.appearance.body{nb+5*(leg-1)} = ...
        {};
    
    % Hip Rx
    nb = 8;
    cyl1 = [0.5*params.leg_rad 0 0;...
        -0.5*params.leg_rad 0 0];
    cyl2 = [0 0 0;...
        side_sign(:,leg)'.*params.hipRyLocation];
    model.appearance.body{nb+5*(leg-1)} = ...
        {};
    
    % Hip Ry
    nb = 9;
    cyl1 = [0 0.5*params.leg_rad 0;...
        0 -0.5*params.leg_rad 0];
    cyl2 = [0 0 0;...
        side_sign(:,leg)'.*params.kneeLocation];
    model.appearance.body{nb+5*(leg-1)} = ...
        {};
    
    % Knee
    nb = 10;
    cyl1 = [0 0.5*params.leg_rad 0;...
        0 -0.5*params.leg_rad 0];
    cyl2 = [0 0 0;...
        side_sign(:,leg)'.*params.ankleLocation];
    model.appearance.body{nb+5*(leg-1)} = ...
        {};
    
    % Ankle
    nb = 11;
    cyl1 = [0 0.5*params.leg_rad 0;...
        0 -0.5*params.leg_rad 0];
    cyl2 = [0 0 0;...
        side_sign(:,leg)'.*[params.footToeLength 0 -params.footHeight]];
    cyl3 = [0 0 0;...
        side_sign(:,leg)'.*[-params.footHeelLength 0 -params.footHeight]];
    model.appearance.body{nb+5*(leg-1)} = ...
        {};
end

for arm = 1:model.NARMS
    % Shoulder Rx
    nb = 17;
    cyl1 = [0.75*params.leg_rad 0 0;...
        -0.75*params.leg_rad 0 0];
    cyl2 = [0 0 0;...
        side_sign(:,arm)'.*params.shoulderRyLocation];
    model.appearance.body{nb+3*(arm-1)} = ...
        {};
    
    % Shoulder Ry
    nb = 18;
    cyl1 = [0 0.75*params.leg_rad 0;...
        0 -0.75*params.leg_rad 0];
    cyl2 = [0 0 0;...
        side_sign(:,arm)'.*params.elbowLocation];
    model.appearance.body{nb+3*(arm-1)} = ...
        {};
    
    % Elbow
    nb = 19;
    cyl1 = [0 0.75*params.leg_rad 0;...
        0 -0.75*params.leg_rad 0];
    cyl2 = [0 0 0;...
        side_sign(:,arm)'.*[0 0 -params.lowerArmLength]];
    model.appearance.body{nb+3*(arm-1)} = ...
        {};
    
end

% Base
model.appearance.base = ...
    { 'box', [-1.5 -1.5 ground_height;...
    10.5 1.5 ground_height-0.01]};
