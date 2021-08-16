%% cleanup
clear all; clc; close all;
restoredefaultpath;

%% setup
addpath(genpath('../../../utilities_general'));
addpath(genpath('../'));
import casadi.*

%% flags
run_IK = false;
show_animation = true;

%% load function
f = Function.load('landingCtrller_KNITRO.casadi');

%% parameters
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModel(params, model);

sideSign = [1 -1 1, 1 1 1, -1 -1 1, -1 1 1];

q_leg_home = [0 -1.45 2.65];
q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
[~,Ibody_val] = get_mass_matrix(model, q_home, 0);
mass_val = Ibody_val(6,6);
Ibody_inv_val = inv(Ibody_val(1:3,1:3));

N = 21; % N = 11
T = 0.6; % T = 0.22
dt_val = repmat(T/(N-1),1,N-1);

q_init_val = [0 0 0 (.25)*(2*rand(1)-1) (pi/4)*(2*rand(1)-1) (.25)*(2*rand(1)-1)]';
qd_init_val = [0.5*(2*rand(1,3)-1) 1*(2*rand(1, 2)-1) -3.5*rand(1)-1]';


for leg = 1:4
    hip_world(:, leg) = rpyToRotMat_xyz(q_init_val(4:6))*params.hipSrbmLocation(leg, :)';
end
td_hip_z = abs(min(hip_world(3,:)));

td_nom = 0.35;
z_max_td = td_nom + td_hip_z + abs(dt_val(1)*qd_init_val(6));

q_init_val(3) = z_max_td;

q_min_val = [-10 -10 0.10 -10 -10 -10];
q_max_val = [10 10 1.0 10 10 10];
qd_min_val = [-10 -10 -10 -40 -40 -40];
qd_max_val = [10 10 10 40 40 40];

q_term_min_val = [-10 -10 0.15 -0.1 -0.1 -10];
q_term_max_val = [10 10 5 0.1 0.1 10];
qd_term_min_val = [-10 -10 -10 -40 -40 -40];
qd_term_max_val = [10 10 10 40 40 40];

q_term_ref = [0 0 0.25, 0 0 0]';
qd_term_ref = [0 0 0, 0 0 0]';

c_ref = diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -0.35],1,4)';
f_ref = zeros(12,1);

QX_val = [0 0 0, 10 10 0, 10 10 10, 10 10 10]';
QX_val = zeros(12, 1);
QN_val = [0 0 100, 100 100 0, 10 10 10, 10 10 10]';
Qc_val = [0 0 0]';
Qf_val = [0.0001 0.0001 0.001]';

mu_val = .5;
l_leg_max_val = .40;
f_max_val = 500;

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

c_init_val = zeros(12, 1);
for leg = 1:4
    xyz_idx = 3*leg-2 : 3*leg;
    p_foot_rel = sideSign(xyz_idx)'.*[0.2 0.15 -0.3]';
    c_init_val(xyz_idx) = q_init_val(1:3) + rpyToRotMat_xyz(q_init_val(4:6))*p_foot_rel;
end

jpos_min_val = repmat([-pi/3, -pi/2, 0]', 4, 1);
jpos_max_val = repmat([pi/3, pi/2, 3*pi/4]', 4, 1);
jpos_guess = repmat([0, -pi/4, pi/2]', 4*(N-1), 1);
v_body = rpyToRotMat_xyz(q_init_val(4:6))'*(qd_init_val(4:6));

kin_box_val = [kin_box_limits(v_body(1), 'x'); kin_box_limits(v_body(2), 'y')];



%% solve optimization

tic
% solve problem by calling f with numerial arguments (for verification)
disp_box('Solving Problem with Solver, c code and simple bounds');
[res.x,res.f] = f(Xref_val, Uref_val,...
        dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
        q_init_val, qd_init_val, c_init_val, ...
        q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
        QN_val, [Xref_val(:);Uref_val(:); jpos_guess],...
        jpos_min_val, jpos_max_val, kin_box_val, mu_val, l_leg_max_val, mass_val,...
        diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
toc
    
% Decompose solution
X_tmp = zeros(12, N);
jpos_tmp = zeros(3*model.NLEGS, N-1);
U_tmp = zeros(6*model.NLEGS, N-1);

res.x = full(res.x);
X_star = reshape(res.x(1:numel(X_tmp)),size(X_tmp));
jpos_star = reshape(res.x(numel(X_tmp) + 1:numel(X_tmp) + numel(jpos_tmp)), size(jpos_tmp));
U_star = reshape(res.x(numel(X_tmp) + numel(jpos_tmp) + 1:numel(X_tmp) + numel(jpos_tmp) + numel(U_tmp)), size(U_tmp));

q_star(1:6,:) = X_star(1:6,:); qd_star(1:6,:) = X_star(7:12,:);
q_star(7:18,1:end-1) = jpos_star; q_star(7:18, end) = q_star(7:18, end-1);

f_star = U_star(13:24, :); p_star = U_star(1:12, :);

t_star = zeros(1,N);
for k = 2:N
    t_star(k) = t_star(k-1) + dt_val(1,k-1);
end

showmotion(model,t_star,q_star)

%% actuator data
% jacobian torque calculation
J_foot = cell(4, N-1);
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
    joint_vel(i, 1:N-2) = diff(jpos_star(i, :))./dt_val(1);
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

make_plots = true;

if make_plots
    t = tiledlayout(4, 3);
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
    plot(t_star(1:end-1), [torque(1, :);torque(4, :);torque(7, :);torque(10, :)], 'r-')
    plot(t_star(1:end-1), model.tauMax(1)*ones(1, N-1), 'r--')
    plot(t_star(1:end-1), -model.tauMax(1)*ones(1, N-1), 'r--')
    plot(t_star(1:end-1), [torque(2, :);torque(5, :);torque(8, :);torque(11, :)], 'g-')
    plot(t_star(1:end-1), model.tauMax(2)*ones(1, N-1), 'g--')
    plot(t_star(1:end-1), -model.tauMax(2)*ones(1, N-1), 'g--')
    plot(t_star(1:end-1), [torque(3, :);torque(6, :);torque(9, :);torque(12, :)], 'b-')
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
        plot(t_star(:), v(i, :))
    end
    xlabel('Time (s)'); ylabel('Voltage (V)')
    axis([0, t_star(end), -26, 26])
    title('Voltage Limits')
    hold off;
    
    
    
%     % Foot locations
%     figure; hold on;
%     for leg = 1:4
%         xyz_idx = 3*leg-2 : 3*leg;
%         for i = 1:N-1
%             R_world_to_body = rpyToRotMat_xyz(q_star(4:6, i))';
%             p_foot_rel(xyz_idx, i) = R_world_to_body*(p_star(xyz_idx, i) - q_star(1:3, i));
%             p_foot_rel(xyz_idx, i) = p_foot_rel(xyz_idx, i) - params.hipSrbmLocation(leg, :)';
%         end
%     end
%     plot(p_foot_rel(1, :), p_foot_rel(2, :))
%     plot(p_foot_rel(4, :), p_foot_rel(5, :))
%     plot(p_foot_rel(7, :), p_foot_rel(8, :))
%     plot(p_foot_rel(10, :), p_foot_rel(11, :))
%     
%     xlabel('x'); ylabel('y'); zlabel('z')
%     title('Foot locations')
%     hold off;
    
end
