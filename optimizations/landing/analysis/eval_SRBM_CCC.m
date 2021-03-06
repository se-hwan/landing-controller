function [X_star, q_star, f_star, p_star] = eval_SRBM_CCC(q_init, qd_init)

    addpath(genpath('../../../utilities_general'));
    addpath(genpath('codegen_casadi'));
    import casadi.*

    run_IK = true;
    
    f = Function.load('codegen_casadi/f_quad_SRBM.casadi');

    disp_box('Building Robot Model');
    params = get_robot_params('mc3D');
    model  = get_robot_model(params);
    model  = buildShowMotionModelMC3D(params, model, 0);

    q_leg_home = [0 -1.45 2.65];
    q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
    [~,Ibody_val] = get_mass_matrix(model, q_home, 0);
    mass_val = Ibody_val(6,6);
    Ibody_inv_val = inv(Ibody_val(1:3,1:3));

    N = 41;
    T = 0.6;
    dt_val = repmat(T/(N-1),1,N-1);
    
    q_init_val = q_init;
    qd_init_val = qd_init;
        
    q_min_val = [-10 -10 0.15 -10 -10 -10];
    q_max_val = [10 10 1.0 10 10 10];
    qd_min_val = [-10 -10 -10 -40 -40 -40];
    qd_max_val = [10 10 10 40 40 40];

    q_term_min_val = [-10 -10 0.15 -0.1 -0.1 -10];
    q_term_max_val = [10 10 5 0.1 0.1 10];
    qd_term_min_val = [-10 -10 -10 -40 -40 -40];
    qd_term_max_val = [10 10 10 40 40 40];

    q_term_ref = [0 0 0.2, 0 0 0]';
    qd_term_ref = [0 0 0, 0 0 0]';

    c_init_val = repmat(q_init_val(1:3),4,1)+...
        diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -q_init_val(3)],1,4)';

    c_ref = diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -0.35],1,4)';
    f_ref = zeros(12,1);

    QX_val = [0 0 0, 10 10 0, 10 10 10, 10 10 10]';
    QX_val = zeros(12, 1);
    QN_val = [0 0 100, 100 100 0, 10 10 10, 10 10 10]';
    Qc_val = [0 0 0]';
    Qf_val = [0.0001 0.0001 0.001]';

    mu_val = 1;
    l_leg_max_val = .35;
    f_max_val = 250;

    for i = 1:6
        Xref_val(i,:)   = linspace(q_init_val(i),q_term_ref(i),N);
        Xref_val(6+i,:) = linspace(qd_init_val(i),qd_term_ref(i),N);
    end
    for leg = 1:4
        for xyz = 1:3
            Uref_val(3*(leg-1)+xyz,:)    = Xref_val(xyz,1:end-1) + c_ref(3*(leg-1)+xyz);
            Uref_val(12+3*(leg-1)+xyz,:) = f_ref(xyz).*ones(1,N-1);
        end
    end

    % evaluate optimization problem with given parameters
    tic
    disp_box('Solving Problem with Solver, c code and simple bounds');
    [res.x,res.f] = f(Xref_val, Uref_val,...
        dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
        q_init_val, qd_init_val, c_init_val,...
        q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
        QX_val, QN_val, Qc_val, Qf_val, [Xref_val(:);Uref_val(:)],...
        mu_val, l_leg_max_val, f_max_val, mass_val,...
        diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
    toc
    
    % Decompose solution
    X = zeros(12, 41);
    U = zeros(6*model.NLEGS, 40);

    res.x = full(res.x);
    X_star = reshape(res.x(1:numel(X)),size(X));
    q_star(1:6,:) = X_star(1:6,:);
    q_star(7:18,:) = repmat(q_home(7:end),1,N);
    U_star = reshape(res.x(numel(X)+1:numel(X)+numel(U)), size(U));
    t_star = zeros(1,N);

    p_star = U_star(1:12, :);
    f_star = U_star(13:24, :);

    for k = 2:N
        t_star(k) = t_star(k-1) + dt_val(1,k-1);
    end
    % inverse kinematics, if called

    q_foot_guess = repmat([0 -0.7 1.45]', 4, 1);
    if run_IK
        for i = 1:N-1
            [x, fval, exitflag] = inverse_kinematics(U_star(1:12,i), model, q_star(1:6,i), q_foot_guess);
            if exitflag <= 0
                q_star(7:18,i) = q_foot_guess;
            end
            q_star(7:18,i) = x;
        end
        q_star(7:18,N) = q_star(7:18,N-1);
    else
        q_star(7:18,:) = repmat(repmat(q_leg_home', 4, 1),1,N);
    end
    
end