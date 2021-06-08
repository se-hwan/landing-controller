function [ model ] = get_robot_leg_model( params )
% Get spatial_v2 model of robot leg

I_l1   = mcI(params.l1_mass, [0 0 -params.l1/2], boxInertia(params.l1_mass,[params.leg_rad*2 params.leg_rad*2 params.l1]));
I_l2   = mcI(params.l2_mass, [0 0 -params.l2/2], boxInertia(params.l2_mass,[params.leg_rad*2 params.leg_rad*2 params.l2]));

NLEGS = 1;
model.NB = 2;
model.NLEGS = 1;
model.gravity = [0 0 -9.81];
model.parent  = zeros(1,model.NB);             % parent body indices
model.jtype   = repmat({'  '},model.NB,1);     % joint types
model.Xtree   = repmat({eye(6)},model.NB,1);   % coordinate transforms
model.I       = repmat({zeros(6)},model.NB,1); % spatial inertias
model.Xfoot   = repmat({zeros(6)},NLEGS,1);    % feet
model.b_foot  = zeros(1,NLEGS);                % body foot is connected to

nb = 0;
%% Hip
nb = nb + 1;
model.parent(nb) = nb - 1;  % hip parent is fixed
model.jtype{nb}  = 'Ry';     % rotate around y
%                 flip so y+ is leg forward,   translate to front/back
model.Xtree{nb}  = plux(rz(pi),[0 0 0]') * plux(eye(3),[params.body_length 0 0]');
model.I{nb}      = I_l1;
    
%% Knee
nb = nb + 1;
model.parent(nb) = nb - 1; % knee parent is hip
model.jtype{nb}  = 'Ry'; % rotate around y
model.Xtree{nb}  = plux(eye(3),[0 0 -params.l1]');
model.I{nb}      = I_l2;

%% Foot (bonus)
model.Xfoot{1}  = plux(eye(3),[0 0 -params.l2]'); 
model.b_foot(1) = nb;

end

% from Pat's model of cheetah 3
function I = boxInertia(mass, x)
    I = (norm(x)^2*eye(3) - diag(x.^2))*mass/12;
end
