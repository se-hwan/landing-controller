clear all; clc; close all;

%% add library paths
addpath(genpath('../../../utilities_general'));
addpath(genpath('../utilties_landing'));
show_plots = true;

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModel(params, model);

%% Load and parse data
input = readtable('data/time_idx_normalization/input_0.csv'); input = input.Var1;
nn_pred = readtable('data/time_idx_normalization/nnpred_0.csv'); nn_pred = nn_pred.Var1;
output = readtable('data/time_idx_normalization/original_0.csv'); output = output.Var1;

N = 21;

[X_nlp, U_nlp, jpos_nlp] = data_denormalization(output);
q_nlp(1:6,:) = X_nlp(1:6, :); qd_nlp = X_nlp(7:12, :);
q_nlp(7:18,1:end-1) = jpos_nlp; q_nlp(7:18, end) = q_nlp(7:18, end-1);
f_nlp = U_nlp(13:24, :); p_nlp = U_nlp(1:12, :);

[X_nn, U_nn, jpos_nn] = data_denormalization(nn_pred);
q_nn(1:6,:) = X_nn(1:6, :); qd_nn = X_nn(7:12, :);
q_nn(7:18,1:end-1) = jpos_nn; q_nn(7:18, end) = q_nn(7:18, end-1);
f_nn = U_nn(13:24, :); p_nn = U_nn(1:12, :);

dt_val = [0.05 repmat(0.02, 1, 15) [0.05 0.05 0.1 0.2]];
t_star = zeros(1,N);
for k = 2:N
    t_star(k) = t_star(k-1) + dt_val(1,k-1);
end

%% visualization - NLP
showmotion_floatingBase(model,t_star(1:end-1),q_nlp(:,1:end-1))
if show_plots
    plot_results(model, params, t_star, X_nlp, U_nlp, jpos_nlp);
end
%% visualization - NN
showmotion_floatingBase(model,t_star(1:end-1),q_nn(:,1:end-1))
if show_plots
    plot_results(model, params, t_star, X_nn, U_nn, jpos_nn);
end

    