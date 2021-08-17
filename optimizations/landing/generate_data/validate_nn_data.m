clear all; clc;

%% add library paths
addpath(genpath('../../../utilities_general'));
make_plots = true;

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModel(params, model);

%% Load and parse data
input = readtable('sehwan_csv/input_0.csv'); input = input.Var1;
nn_pred = readtable('sehwan_csv/nnpred_0.csv'); nn_pred = nn_pred.Var1;
output = readtable('sehwan_csv/output_0.csv'); output = output.Var1;

N = 21;
X_tmp = zeros(12, N);
U_tmp = zeros(6*model.NLEGS, N-1);
jpos_tmp = zeros(3*model.NLEGS, N-1);

X_nlp = reshape(output(1:numel(X_tmp)),size(X_tmp));
U_nlp = reshape(output(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
jpos_nlp = reshape(output(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(U_tmp) + numel(jpos_tmp)), size(jpos_tmp));
q_nlp(1:6,:) = X_nlp(1:6, :); qd_nlp = X_nlp(7:12, :);
q_nlp(7:18,1:end-1) = jpos_nlp; q_nlp(7:18, end) = q_nlp(7:18, end-1);
f_nlp = U_nlp(13:24, :); p_nlp = U_nlp(1:12, :);

X_nn = reshape(nn_pred(1:numel(X_tmp)),size(X_tmp));
U_nn = reshape(nn_pred(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
jpos_nn = reshape(nn_pred(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(U_tmp) + numel(jpos_tmp)), size(jpos_tmp));
q_nn(1:6,:) = X_nn(1:6, :); qd_nn = X_nn(7:12, :);
q_nn(7:18,1:end-1) = jpos_nn; q_nn(7:18, end) = q_nn(7:18, end-1);
f_nn = U_nlp(13:24, :); p_nn = U_nlp(1:12, :);

t_star = zeros(1,N);
for k = 2:N
    t_star(k) = t_star(k-1) + 0.03;
end

%% visualization - NLP
showmotion(model,t_star,q_nlp)

%% visualization - NN
showmotion(model,t_star,q_nn)

%% actuator data
% jacobian torque calculation
torque_nlp = zeros(12, N-1);
for i = 1:N-1
    R_world_to_body = rpyToRotMat_xyz(q_nlp(4:6, i))';
    J_f = get_foot_jacobians_mc(model, params, jpos_nlp(:, i));
    for leg = 1:4
        xyz_idx = 3*leg-2:3*leg;
        torque_nlp(xyz_idx, i) = J_f{leg}'*(-R_world_to_body*f_nlp(xyz_idx, i));
    end
end

% motor voltage calculation
v_nlp = zeros(12, N-1);
joint_vel_nlp = zeros(12, N-1);
for i = 1:12
    joint_vel_nlp(i, 1:N-2) = diff(jpos_nlp(i, :))./0.03;
end

for i = 1:N-1
    tau_motor_des_i = torque_nlp(:,i) ./ repmat(model.gr,4,1);
    current_des_i = tau_motor_des_i ./ (1.5*repmat(model.kt, 4, 1));
    joint_vel_i = joint_vel_nlp(:, i);
    back_emf_i = joint_vel_i .* repmat(model.gr, 4, 1) .* repmat(model.kt, 4, 1) * 2.0;
    v_des_i = current_des_i .* repmat(model.Rm, 4, 1) + back_emf_i;
    v_nlp(:, i) = v_des_i;
end
v_nlp = [zeros(12, 1) v_nlp];

torque_nn = zeros(12, N-1);
for i = 1:N-1
    R_world_to_body = rpyToRotMat_xyz(q_nn(4:6, i))';
    J_f = get_foot_jacobians_mc(model, params, jpos_nn(:, i));
    for leg = 1:4
        xyz_idx = 3*leg-2:3*leg;
        torque_nn(xyz_idx, i) = J_f{leg}'*(-R_world_to_body*f_nn(xyz_idx, i));
    end
end

% motor voltage calculation
v_nn = zeros(12, N-1);
joint_vel_nn = zeros(12, N-1);
for i = 1:12
    joint_vel_nn(i, 1:N-2) = diff(jpos_nn(i, :))./0.03;
end

for i = 1:N-1
    tau_motor_des_i = torque_nn(:,i) ./ repmat(model.gr,4,1);
    current_des_i = tau_motor_des_i ./ (1.5*repmat(model.kt, 4, 1));
    joint_vel_i = joint_vel_nn(:, i);
    back_emf_i = joint_vel_i .* repmat(model.gr, 4, 1) .* repmat(model.kt, 4, 1) * 2.0;
    v_des_i = current_des_i .* repmat(model.Rm, 4, 1) + back_emf_i;
    v_nn(:, i) = v_des_i;
end
v_nn = [zeros(12, 1) v_nn];


%% plots
if make_plots
    figure; t = tiledlayout(4, 3);
    t.Padding = 'compact';
    t.TileSpacing = 'compact';
    
    % GRFs
    nexttile
    hold on;
    for leg = 1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        plot(t_star(1:end-1), f_nlp(xyz_idx(3), :));
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
        plot(t_star(1:end-1), f_nlp(xyz_idx(1), :));
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
        plot(t_star(1:end-1), f_nlp(xyz_idx(2), :));
    end
    xlabel('Time (s)'); ylabel('Force (N)');
    title('Y ground reaction forces');
    legend('FR', 'FL', 'BR', 'BL')
    hold off;
    
    % CoM posn
    nexttile
    hold on;
    plot(t_star, q_nlp(1,:))
    plot(t_star, q_nlp(2,:))
    plot(t_star, q_nlp(3,:))
    xlabel('Time (s)'); ylabel('Position (m)');
    legend('X','Y','Z')
    title('CoM Position')
    hold off;
    
    % CoM posn
    nexttile
    hold on;
    plot(t_star, qd_nlp(4,:))
    plot(t_star, qd_nlp(5,:))
    plot(t_star, qd_nlp(6,:))
    xlabel('Time (s)'); ylabel('Velocity (m/s)');
    legend('X','Y','Z')
    title('CoM Velocity')
    hold off;
    
    % orientation
    nexttile
    hold on;
    plot(t_star, rad2deg(q_nlp(4,:)))
    plot(t_star, rad2deg(q_nlp(5,:)))
    plot(t_star, rad2deg(q_nlp(6,:)))
    xlabel('Time (s)'); ylabel('Orientation (degrees)');
    legend('Roll', 'Pitch', 'Yaw')
    title('CoM Orientation')
    hold off;
    
    % torque limits
    nexttile([1 3])
    hold on;
    plot(t_star(1:end-1), [torque_nlp(1, :);torque_nlp(4, :);torque_nlp(7, :);torque_nlp(10, :)], 'r-')
    plot(t_star(1:end-1), model.tauMax(1)*ones(1, N-1), 'r--')
    plot(t_star(1:end-1), -model.tauMax(1)*ones(1, N-1), 'r--')
    plot(t_star(1:end-1), [torque_nlp(2, :);torque_nlp(5, :);torque_nlp(8, :);torque_nlp(11, :)], 'g-')
    plot(t_star(1:end-1), model.tauMax(2)*ones(1, N-1), 'g--')
    plot(t_star(1:end-1), -model.tauMax(2)*ones(1, N-1), 'g--')
    plot(t_star(1:end-1), [torque_nlp(3, :);torque_nlp(6, :);torque_nlp(9, :);torque_nlp(12, :)], 'b-')
    plot(t_star(1:end-1), model.tauMax(3)*ones(1, N-1), 'b--')
    plot(t_star(1:end-1), -model.tauMax(3)*ones(1, N-1), 'b--')
    xlabel('Time (s)'); ylabel('Torque (Nm)')
    title("Torque Limits")
    hold off;
    
        
    % voltage limits
    nexttile([1 3])
    hold on;
    plot(t_star(:), model.batteryV*ones(1, N), 'k--')
    plot(t_star(:), -model.batteryV*ones(1, N), 'k--')
    for i = 1:12
        plot(t_star(:), v_nlp(i, :))
    end
    xlabel('Time (s)'); ylabel('Voltage (V)')
    axis([0, t_star(end), -26, 26])
    title('Voltage Limits')
    hold off;
    
    
    
    figure; t_nn = tiledlayout(4, 3);
    t_nn.Padding = 'compact';
    t_nn.TileSpacing = 'compact';
    
    % GRFs
    nexttile
    hold on;
    for leg = 1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        plot(t_star(1:end-1), f_nn(xyz_idx(3), :));
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
        plot(t_star(1:end-1), f_nn(xyz_idx(1), :));
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
        plot(t_star(1:end-1), f_nn(xyz_idx(2), :));
    end
    xlabel('Time (s)'); ylabel('Force (N)');
    title('Y ground reaction forces');
    legend('FR', 'FL', 'BR', 'BL')
    hold off;
    
    % CoM posn
    nexttile
    hold on;
    plot(t_star, q_nn(1,:))
    plot(t_star, q_nn(2,:))
    plot(t_star, q_nn(3,:))
    xlabel('Time (s)'); ylabel('Position (m)');
    legend('X','Y','Z')
    title('CoM Position')
    hold off;
    
    % CoM posn
    nexttile
    hold on;
    plot(t_star, qd_nn(4,:))
    plot(t_star, qd_nn(5,:))
    plot(t_star, qd_nn(6,:))
    xlabel('Time (s)'); ylabel('Velocity (m/s)');
    legend('X','Y','Z')
    title('CoM Velocity')
    hold off;
    
    % orientation
    nexttile
    hold on;
    plot(t_star, rad2deg(q_nn(4,:)))
    plot(t_star, rad2deg(q_nn(5,:)))
    plot(t_star, rad2deg(q_nn(6,:)))
    xlabel('Time (s)'); ylabel('Orientation (degrees)');
    legend('Roll', 'Pitch', 'Yaw')
    title('CoM Orientation')
    hold off;
    
    % torque limits
    nexttile([1 3])
    hold on;
    plot(t_star(1:end-1), [torque_nn(1, :);torque_nn(4, :);torque_nn(7, :);torque_nn(10, :)], 'r-')
    plot(t_star(1:end-1), model.tauMax(1)*ones(1, N-1), 'r--')
    plot(t_star(1:end-1), -model.tauMax(1)*ones(1, N-1), 'r--')
    plot(t_star(1:end-1), [torque_nn(2, :);torque_nn(5, :);torque_nn(8, :);torque_nn(11, :)], 'g-')
    plot(t_star(1:end-1), model.tauMax(2)*ones(1, N-1), 'g--')
    plot(t_star(1:end-1), -model.tauMax(2)*ones(1, N-1), 'g--')
    plot(t_star(1:end-1), [torque_nn(3, :);torque_nn(6, :);torque_nn(9, :);torque_nn(12, :)], 'b-')
    plot(t_star(1:end-1), model.tauMax(3)*ones(1, N-1), 'b--')
    plot(t_star(1:end-1), -model.tauMax(3)*ones(1, N-1), 'b--')
    xlabel('Time (s)'); ylabel('Torque (Nm)')
    title("Torque Limits")
    hold off;
    
        
    % voltage limits
    nexttile([1 3])
    hold on;
    plot(t_star(:), model.batteryV*ones(1, N), 'k--')
    plot(t_star(:), -model.batteryV*ones(1, N), 'k--')
    for i = 1:12
        plot(t_star(:), v_nn(i, :))
    end
    xlabel('Time (s)'); ylabel('Voltage (V)')
    axis([0, t_star(end), -26, 26])
    title('Voltage Limits')
    hold off;
end

    
    