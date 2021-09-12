clear; clc; close all;

%% add library paths
addpath(genpath('../../../utilities_general'));
addpath(genpath('../utilities_landing'));
import casadi.*
show_plots = false;
show_animation = false;

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModel(params, model);

%% import trained network
modelfile = 'neural_networks/nn_TO_landing.onnx';
nn_params = importONNXFunction(modelfile,'nn_TO_landing');

modelfile = 'neural_networks/old/indexedTime_2x64_bs16.onnx';
nn_params = importONNXFunction(modelfile,'indexedTime_2x64_bs16');

modelfile = 'neural_networks/nn-landing-NODES-64-BATCHSIZE-16.onnx';
nn_params = importONNXFunction(modelfile,'nn_landing_n64_b16');

%% load normalized statistics
load('data/landing_norm_param.mat');
input_mean = data_stats.mean_input;
input_std = data_stats.std_input;

%% load functions
f_ipopt_SRBM = Function.load('../codegen_casadi/landingCtrller_IPOPT.casadi');
f_knitro = Function.load('../codegen_casadi/landingCtrller_KNITRO.casadi');
f_knitro_ws = Function.load('../codegen_casadi/landingCtrller_KNITRO_ws.casadi');
% f_knitro_ws_MA97 = Function.load('../codegen_casadi/landingCtrller_KNITRO_ws_MA97.casadi');

%% begin trials

N_trials = 1;
t_solve = zeros(4, N_trials); % rows: nn ws eval time, nn ws solve time, KD cold start, SRBM warm start

N = 21; 
X_tmp = zeros(12, N);
U_tmp = zeros(6*model.NLEGS, N-1);
jpos_tmp = zeros(3*model.NLEGS, N-1);

for n = 1:N_trials
    
%% initial falling conditions
q_init_val = [0 0 0 (0.25)*(2*rand(1)-1) (0.25)*(2*rand(1)-1) (.25)*(2*rand(1)-1)]';     % major pitch
qd_init_val = [0.1*(2*rand(1,3)-1) 1*(2*rand(1, 2)-1) -3*rand(1)-3]';

q_init_val = [0 0 0 (0.25)*(2*rand(1)-1) pi/4 (.25)*(2*rand(1)-1)]';     % major pitch
qd_init_val = [0.1*(2*rand(1,3)-1) 1*(2*rand(1, 2)-1) -3*rand(1)-3]';

q_init_val = [0 0 0 0 pi/4 pi/4]'; qd_init_val = [0 0 0 0 0 -3]';
%% timestep parameters
dt_val = [0.05 repmat(0.02, 1, 15) [0.05 0.05 0.1 0.2]];
t_star = zeros(1,N);
for k = 2:N
    t_star(k) = t_star(k-1) + dt_val(1,k-1);
end

%% parameters
sideSign = [1 -1 1, 1 1 1, -1 -1 1, -1 1 1];

for leg = 1:4
    hip_world(:, leg) = rpyToRotMat(q_init_val(4:6))*params.hipSrbmLocation(leg, :)';
end
td_hip_z = abs(min(hip_world(3,:)));

td_nom = 0.35;
z_max_td = td_nom + td_hip_z + abs(dt_val(1)*qd_init_val(6));

q_init_val(3) = z_max_td;

q_term_min_val = [-10 -10 0.15 -0.1 -0.1 -10];
q_term_max_val = [10 10 5 0.1 0.1 10];
qd_term_min_val = [-10 -10 -10 -.5 -.5 -.5];
qd_term_max_val = [10 10 10 .5 .5 .5];

q_min_val = [-10 -10 0.075 -10 -10 -10];
q_max_val = [10 10 1.0 10 10 10];
qd_min_val = [-10 -10 -10 -40 -40 -40];
qd_max_val = [10 10 10 40 40 40];

q_term_ref = [0 0 0.225, 0 0 0]';
qd_term_ref = [0 0 0, 0 0 0]';

c_init_val = zeros(12, 1);
for leg = 1:4
    xyz_idx = 3*leg-2 : 3*leg;
    p_foot_rel = sideSign(xyz_idx)'.*[0.2 0.15 -0.3]';
    c_init_val(xyz_idx) = q_init_val(1:3) + rpyToRotMat(q_init_val(4:6))*p_foot_rel;
end

q_leg_home = [0 -1.45 2.65];
q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
[~,Ibody_val] = get_mass_matrix(model, q_home, 0);
mass_val = Ibody_val(6,6);
Ibody_inv_val = inv(Ibody_val(1:3,1:3));

jpos_min_val = repmat([-pi/3, -pi/2, 0]', 4, 1);
jpos_max_val = repmat([pi/3, pi/2, 3*pi/4]', 4, 1);

v_body = rpyToRotMat(q_init_val(4:6))'*(qd_init_val(4:6));

kin_box_val = [kin_box_limits(v_body(1), 'x'); kin_box_limits(v_body(2), 'y')];

c_ref = diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.15 -0.3],1,4)';
f_ref = zeros(12,1);

QN_val = [0 0 100, 10 10 0, 10 10 10, 10 10 10]';

mu_val = 0.75;
l_leg_max_val = .4;
f_max_val = 500;

for i = 1:6
    Xref_val(i,:)   = linspace(q_init_val(i),q_term_ref(i),N);
    Xref_val(6+i,:) = linspace(qd_init_val(i),qd_term_ref(i),N);
end
for leg = 1:4
    for xyz = 1:3
        Uref_val(12+3*(leg-1)+xyz,:) = f_ref(xyz).*ones(1,N-1);
    end
end
for i = 1:N-1
    for leg = 1:4
        xyz_idx = 3*leg-2 : 3*leg;
        Uref_val(xyz_idx, i) = Xref_val(1:3, i) + rpyToRotMat(Xref_val(4:6, i))*c_ref(xyz_idx);
    end
end

X_star_guess = Xref_val(:);
U_star_guess = Uref_val(:);
jpos_star_guess = repmat([0, -pi/4, pi/2]', 4*(N-1), 1);

%% load and run neural network

nn_input = [q_init_val(4:6)' qd_init_val'];
nn_input = (nn_input - input_mean')./input_std';

inputPerm = preparePermutationVector(["FeaturesLength","SequenceLength","BatchSize"],...
    ["SequenceLength","BatchSize","FeaturesLength"]);

tic
% nn_ws = double(nn_landing_n64_b16(nn_input', nn_params, 'InputDataPermutation', inputPerm));
nn_ws = double(indexedTime_2x64_bs16(nn_input', nn_params, 'InputDataPermutation', inputPerm));
% nn_ws = double(indexedTime_2x64_bs16(nn_input, nn_params));

[X_nn, U_nn, jpos_nn] = data_denormalization(nn_ws);
% U_nn(13:24, :) = zeros(12, N-1);
q_nn(1:6,:) = X_nn(1:6, :); qd_nn = X_nn(7:12, :);
q_nn(7:18,1:end-1) = jpos_nn; q_nn(7:18, end) = q_nn(7:18, end-1);
f_nn = U_nn(13:24, :); p_nn = U_nn(1:12, :);

if show_animation
    showmotion_floatingBase(model,t_star(1:end-1),q_nn(:,1:end-1))
    plot_results(model, params, t_star, X_nn, U_nn, jpos_nn);
end

%% solve with NN warm start
% tic
[res.x,res.f] = f_knitro_ws(Xref_val, Uref_val,...
        dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
        q_init_val, qd_init_val, c_init_val, ...
        q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
        QN_val, [X_nn(:); U_nn(:); jpos_nn(:)],...
        jpos_min_val, jpos_max_val, kin_box_val, mu_val, l_leg_max_val, mass_val,...
        diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
t_solve_nn_ws = toc; disp("NN WARM STARTED SOLVE TIME: " + num2str(t_solve_nn_ws));
t_solve(2,n) = t_solve_nn_ws;

% U_nn(13:24, :) = zeros(12, N-1);
[res.x,res.f] = f_knitro(Xref_val, Uref_val,...
        dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
        q_init_val, qd_init_val, c_init_val, ...
        q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
        QN_val, [X_nn(:); U_nn(:); jpos_nn(:)],...
        jpos_min_val, jpos_max_val, kin_box_val, mu_val, l_leg_max_val, mass_val,...
        diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
t_solve_nn_ws = toc; disp("NN WARM STARTED SOLVE TIME v2: " + num2str(t_solve_nn_ws));
t_solve(2,n) = t_solve_nn_ws;


[res.x,res.f] = f_knitro_ws(Xref_val, Uref_val,...
        dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
        q_init_val, qd_init_val, c_init_val, ...
        q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
        QN_val, [X_nn(:); U_nn(:); jpos_nn(:)],...
        jpos_min_val, jpos_max_val, kin_box_val, mu_val, l_leg_max_val, mass_val,...
        diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
t_solve_nn_ws = toc; disp("NN WARM STARTED SOLVE TIME ZERO FORCES: " + num2str(t_solve_nn_ws));
t_solve(2,n) = t_solve_nn_ws;

res.x = full(res.x);
X_star = reshape(res.x(1:numel(X_tmp)),size(X_tmp));
U_star = reshape(res.x(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
jpos_star = reshape(res.x(numel(X_tmp) + + numel(U_tmp) + 1:numel(X_tmp) + numel(U_tmp) + numel(jpos_tmp)), size(jpos_tmp));

q_star(1:6,:) = X_star(1:6,:); qd_star(1:6,:) = X_star(7:12,:);
q_star(7:18,1:end-1) = jpos_star; q_star(7:18, end) = q_star(7:18, end-1);

f_star = U_star(13:24, :); p_star = U_star(1:12, :);


if show_animation
    showmotion_floatingBase(model,t_star(1:end-1),q_star(:,1:end-1))
end
if show_plots
    plot_results(model, params, t_star, X_star, U_star, jpos_star);
end


%% cold start
tic
[res.x,res.f] = f_knitro(Xref_val, Uref_val,...
        dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
        q_init_val, qd_init_val, c_init_val, ...
        q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
        QN_val, [X_star_guess(:);jpos_star_guess; U_star_guess(:)],...
        jpos_min_val, jpos_max_val, kin_box_val, mu_val, l_leg_max_val, mass_val,...
        diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
t_solve_cs = toc; disp("COLD STARTED SOLVE TIME: " + num2str(t_solve_cs))

%% initial guess (SRBM LCP solution)
tic
[res.x,res.f] = f_ipopt_SRBM(Xref_val, Uref_val,...
    dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
    q_init_val, qd_init_val, ...
    q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
    QN_val, [Xref_val(:);Uref_val(:)],...
    mu_val, l_leg_max_val, f_max_val, mass_val,...
    diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));

res.x = full(res.x);
X_star_guess = reshape(res.x(1:numel(X_tmp)),size(X_tmp));
U_star_guess = reshape(res.x(numel(X_tmp)+1:numel(X_tmp)+numel(U_tmp)), size(U_tmp));
jpos_star_guess = repmat([0, -pi/4, pi/2]', 4*(N-1), 1);

[res.x,res.f] = f_knitro_ws(Xref_val, Uref_val,...
        dt_val,q_min_val, q_max_val, qd_min_val, qd_max_val,...
        q_init_val, qd_init_val, c_init_val, ...
        q_term_min_val, q_term_max_val, qd_term_min_val, qd_term_max_val,...
        QN_val, [X_star_guess(:);jpos_star_guess; U_star_guess(:)],...
        jpos_min_val, jpos_max_val, kin_box_val, mu_val, l_leg_max_val, mass_val,...
        diag(Ibody_val(1:3,1:3)), diag(Ibody_inv_val(1:3,1:3)));
t_solve_ws = toc; disp("WARM STARTED SOLVE TIME: " + num2str(t_solve_ws))

disp("Initial position conditions: [" + num2str(q_init_val') + "]");
disp("Initial velocity conditions: [" + num2str(qd_init_val') + "]");

t_solve(:, n) = [0; t_solve_nn_ws; t_solve_cs; t_solve_ws];

end


%% visualization
if show_animation
    showmotion(model,t_star(1:end-1),q_star(:,1:end-1))
end
if show_plots
    plot_results(model, params, t_star, X_star, U_star, jpos_star);
end

%%
figure; hold on;
boxplot([t_solve(2,:)' t_solve(3,:)' t_solve(4,:)'], {"NN WS", "CS", "SRBM WS"});
ylabel('Solve times (s)')
hold off;

