function model = buildShowMotionModelMC3D(params, model)

body_c = [204, 126, 31]./255;
l_gray = [130, 153, 181]./255;
black = [28, 28, 28]./255;
d_gray = [125, 125, 125]./255;

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
    {'colour',body_c,...
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
        {'colour',d_gray,...
        'cyl', cyl1, 1.65*params.leg_rad,...
        'colour',l_gray,...
        'cyl', cyl2, params.leg_rad};
    
    % Hip
    nb = 8;
    cyl1 = [0 0.5*params.leg_rad 0;...
        0 -0.5*params.leg_rad 0];
    cyl2 = [0 0 0;...
        side_sign(:,leg)'.*params.kneeLocation];
    model.appearance.body{nb+3*(leg-1)} = ...
        {'colour',d_gray,...
        'cyl', cyl1, 1.45*params.leg_rad,...
        'colour',l_gray,...
        'cyl', cyl2, params.leg_rad};
    
    % Knee
    nb = 9;
    cyl1 = [0 0.5*params.leg_rad 0;...
        0 -0.5*params.leg_rad 0];
    cyl2 = [0 0 0;...
        side_sign(:,leg)'.*params.footLocation];
    model.appearance.body{nb+3*(leg-1)} = ...
        {'colour',d_gray,...
        'cyl', cyl1, 1.25*params.leg_rad,...
        'colour',l_gray,...
        'cyl', cyl2, params.leg_rad,...
        'colour',black,...
        'sphere',side_sign(:,leg)'.*params.footLocation, 1.05*params.leg_rad}; 
end