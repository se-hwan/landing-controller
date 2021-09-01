function [] = plot_results(model, params, t_star, X_star, U_star, jpos_star)
    q_star = X_star(1:6, :);
    qd_star = X_star(7:12, :);
    p_star = U_star(1:12, :);
    f_star = U_star(13:24, :);
    q_star(7:18,1:end-1) = jpos_star;
    q_star(7:18, end) = q_star(7:18, end-1);
    
    N = length(t_star);
    dt_val = diff(t_star);
    
    %% actuator data
    % jacobian torque calculation
    torque = zeros(12, N-1);
    for i = 1:N-1
        R_world_to_body = rpyToRotMat_xyz(q_star(4:6, i))';
        J_f = get_foot_jacobians_mc(model, params, jpos_star(:, i));
        for leg = 1:4
            xyz_idx = 3*leg-2:3*leg;
            torque(xyz_idx, i) = J_f{leg}'*(-R_world_to_body*f_star(xyz_idx, i));
        end
    end

    % motor voltage calculation
    v = zeros(12, N-1);
    joint_vel = zeros(12, N-1);
    for i = 1:12
        joint_vel(i, 1:N-2) = diff(jpos_star(i, :))./dt_val(1:N-2);
    end

    for i = 1:N-1
        tau_motor_des_i = torque(:,i) ./ repmat(model.gr,4,1);
        current_des_i = tau_motor_des_i ./ (1.5*repmat(model.kt, 4, 1));
        joint_vel_i = joint_vel(:, i);
        back_emf_i = joint_vel_i .* repmat(model.gr, 4, 1) .* repmat(model.kt, 4, 1) * 2.0;
        v_des_i = current_des_i .* repmat(model.Rm, 4, 1) + back_emf_i;
        v(:, i) = v_des_i;
    end
    v = [zeros(12, 1) v];

    %% plots
    figure;
    t = tiledlayout(3, 3);
    t.Padding = 'compact';
    t.TileSpacing = 'compact';

    % GRFs
    nexttile
    hold on;
    for leg = 1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        plot(t_star(1:end-1), f_star(xyz_idx(3), :));
    end
    xlabel('Time (s)'); ylabel('Force (N)');
    title('Vertical ground reaction forces');
    legend('FR', 'FL', 'BR', 'BL')
    hold off;
    
    % X GRFs
    nexttile
    hold on;
    for leg = 1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        plot(t_star(1:end-1), f_star(xyz_idx(1), :));
    end
    xlabel('Time (s)'); ylabel('Force (N)');
    title('X ground reaction forces');
    legend('FR', 'FL', 'BR', 'BL')
    hold off;

    % Y GRFs
    nexttile
    hold on;
    for leg = 1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        plot(t_star(1:end-1), f_star(xyz_idx(2), :));
    end
    xlabel('Time (s)'); ylabel('Force (N)');
    title('Y ground reaction forces');
    legend('FR', 'FL', 'BR', 'BL')
    hold off;

    % CoM posn
    nexttile
    hold on;
    plot(t_star, q_star(1,:))
    plot(t_star, q_star(2,:))
    plot(t_star, q_star(3,:))
    xlabel('Time (s)'); ylabel('Position (m)');
    legend('X','Y','Z')
    title('CoM Position')
    hold off;

    % CoM posn
    nexttile
    hold on;
    plot(t_star, qd_star(4,:))
    plot(t_star, qd_star(5,:))
    plot(t_star, qd_star(6,:))
    xlabel('Time (s)'); ylabel('Velocity (m/s)');
    legend('X','Y','Z')
    title('CoM Velocity')
    hold off;

    % orientation
    nexttile
    hold on;
    plot(t_star, rad2deg(q_star(4,:)))
    plot(t_star, rad2deg(q_star(5,:)))
    plot(t_star, rad2deg(q_star(6,:)))
    xlabel('Time (s)'); ylabel('Orientation (degrees)');
    legend('Roll', 'Pitch', 'Yaw')
    title('CoM Orientation')
    hold off;

    % torque limits
    nexttile([1 3])
    hold on;
    L(1) = plot(0, 0, 'r-');
    plot(t_star(1:end-1), [torque(1, :);torque(4, :);torque(7, :);torque(10, :)], 'r-');
    plot(t_star(1:end-1), model.tauMax(1)*ones(1, N-1), 'r--');
    plot(t_star(1:end-1), -model.tauMax(1)*ones(1, N-1), 'r--');
    L(2) = plot(0, 0, 'g-');
    plot(t_star(1:end-1), [torque(2, :);torque(5, :);torque(8, :);torque(11, :)], 'g-');
    plot(t_star(1:end-1), model.tauMax(2)*ones(1, N-1), 'g--');
    plot(t_star(1:end-1), -model.tauMax(2)*ones(1, N-1), 'g--');
    L(3) = plot(0, 0, 'b-');
    plot(t_star(1:end-1), [torque(3, :);torque(6, :);torque(9, :);torque(12, :)], 'b-');
    plot(t_star(1:end-1), model.tauMax(3)*ones(1, N-1), 'b--');
    plot(t_star(1:end-1), -model.tauMax(3)*ones(1, N-1), 'b--');
    xlabel('Time (s)'); ylabel('Torque (Nm)')
    legend(L, {'Ab/ad', 'Hip', 'Knee'})
    title("Torque Limits")
    hold off;

%     % voltage limits
%     nexttile([1 3])
%     hold on;
%     plot(t_star(:), model.batteryV*ones(1, N), 'k--')
%     plot(t_star(:), -model.batteryV*ones(1, N), 'k--')
%     for i = 1:12
%         plot(t_star(:), v(i, :))
%     end
%     xlabel('Time (s)'); ylabel('Voltage (V)')
%     axis([0, t_star(end), -26, 26])
%     title('Voltage Limits')
%     hold off;

end