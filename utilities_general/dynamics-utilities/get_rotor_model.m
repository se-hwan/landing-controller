function [ rotor_model ] = get_rotor_model( params )

hip_x = [params.body_length -params.body_length]/2;

I_rotor = 1 *mcI(params.rotor_mass, [0 0 0], params.rotor_inertia);
rotor_model.NR = 4;
rotor_model.mu = zeros(1,rotor_model.NR);
rotor_model.gamma = zeros(1,rotor_model.NR);
rotor_model.gr = zeros(1,rotor_model.NR);
rotor_model.I = repmat({zeros(6)},rotor_model.NR,1);
rotor_model.X_mu = repmat({eye(6)},rotor_model.NR,1);

% front hip
rotor_model.gr(1) = params.gr_abad_hip;
rotor_model.mu(1) = 3; % fixed to body (floating base)
rotor_model.gamma(1) = 4; % rotation constraint body (upper link)
rotor_model.I{1} = I_rotor; % default y axis rotor
rotor_model.X_mu{1} = plux(eye(3),[hip_x(1) 0 0]'); % relative to fixed-to body (floating base)


% front knee
rotor_model.gr(2) = params.gr_knee;
rotor_model.mu(2) = 4; % fixed to body (upper link)
rotor_model.gamma(2) = 5; % rotation constraint  body (lower link)
rotor_model.I{2} = I_rotor; % default y axis rotor
rotor_model.X_mu{2} = eye(6); % knee rotor at upper link

% rear hip
rotor_model.gr(3) = params.gr_abad_hip;
rotor_model.mu(3) = 3; % fixed to body (float base)
rotor_model.gamma(3) = 6; % rotation constraint body (upper link)
rotor_model.I{3} = I_rotor; % default y axis rotor
rotor_model.X_mu{3} = plux(eye(3),[hip_x(2) 0 0]');

% rear knee
rotor_model.gr(4) = params.gr_knee;
rotor_model.mu(4) = 6; % fixed to body (upper link)
rotor_model.gamma(4) = 7; % rotation constraint coordinate (lower link)
rotor_model.I{4} = I_rotor; % default y axis rotor
rotor_model.X_mu{4} = eye(6); % knee rotor at upper link origin


end

