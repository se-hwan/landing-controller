function parameters = generate_optimization_parameters()

    %% parameters
    disp_box('Building Robot Model');
    params = get_robot_params('mc3D');
    model  = get_robot_model(params);
    model  = buildShowMotionModel(params, model);

    q_leg_home = [0 -1.45 2.65];
    q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
    [~,Ibody_val] = get_mass_matrix(model, q_home, 0);
    mass_val = Ibody_val(6,6);
    Ibody_inv_val = inv(Ibody_val(1:3,1:3));

    N = 21; % N = 11
    T = 0.6; % T = 0.22
    dt_val = repmat(T/(N-1),1,N-1);

    q_min_val = [-10 -10 0.10 -10 -10 -10];
    q_max_val = [10 10 1.0 10 10 10];
    qd_min_val = [-10 -10 -10 -40 -40 -40];
    qd_max_val = [10 10 10 40 40 40];

    q_term_min_val = [-10 -10 0.15 -0.1 -0.1 -10];
    q_term_max_val = [10 10 5 0.1 0.1 10];
    qd_term_min_val = [-10 -10 -10 -40 -40 -40];
    qd_term_max_val = [10 10 10 40 40 40];


    QN_val = [0 0 100, 100 100 0, 10 10 10, 10 10 10]';

    mu_val = 1;
    l_leg_max_val = .4;
    f_max_val = 225;

    
    
    Ibody_diag = diag(Ibody_val(1:3,1:3));
    Ibody_inv_diag = diag(Ibody_inv_val(1:3,1:3));
    
    save('optimization_parameters.mat',...
         'N','dt_val','q_min_val', 'q_max_val', 'qd_min_val', 'qd_max_val',...
         'q_term_min_val', 'q_term_max_val', 'qd_term_min_val', 'qd_term_max_val',...
         'QN_val', 'mu_val', 'l_leg_max_val', 'f_max_val', 'mass_val',...
         'Ibody_diag', 'Ibody_inv_diag','model','params', 'q_home');
    
    
end