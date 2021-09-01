%% cleanup
clear; clc; restoredefaultpath;

%% add library paths
addpath(genpath('../../../utilities_general'));
addpath(genpath('../utilities_landing'));

%% flags
show_plots = true;
show_animation = true;

%% build robot model
disp_box('Building Robot Model');
params = get_robot_params('mc3D');
model  = get_robot_model(params);
model  = buildShowMotionModel(params, model);

%% load normalized data
load('data/training_data_normalized_idx.mat')

%% parameters
N = 21; 
X_tmp = zeros(12, N);
U_tmp = zeros(6*model.NLEGS, N-1);
jpos_tmp = zeros(3*model.NLEGS, N-1);

q_leg_home = [0 -1.45 2.65];
q_home = [0 0 0 0 0 0 q_leg_home q_leg_home q_leg_home q_leg_home]';
[~,Ibody_val] = get_mass_matrix(model, q_home, 0);
mass_val = Ibody_val(6,6);

dt_val = [0.05 repmat(0.02, 1, 15) [0.05 0.05 0.1 0.2]];
t_star = zeros(1,N);
for k = 2:N
    t_star(k) = t_star(k-1) + dt_val(1,k-1);
end

%% bad training data
idx_bad_samples = [766, 166, 480, 205, 293, 930, 450, 830, 113, 290, 735, 674];

sample_idx = 766;
sample_idx = 166;
sample_idx = 480;
sample_idx = 205; % not bad, but sus
sample_idx = 293;
sample_idx = 930;
sample_idx = 450; % not terrible
sample_idx = 830;
sample_idx = 113;
sample_idx = 290;
sample_idx = 601; % fine
sample_idx = 664; % fine
sample_idx = 735; % not bad, but unusual
sample_idx = 674;


%% trimmed data
training_data_trimmed = struct();
training_data_trimmed.input = training_data_normalized.input;
training_data_trimmed.output = training_data_normalized.output;

for i = 1:length(idx_bad_samples)
    training_data_trimmed.input(:, idx_bad_samples(i)) = [];
    training_data_trimmed.output(:, idx_bad_samples(i)) = [];
end

save('data/training_data_normalized_idx_trimmed.mat', 'training_data_trimmed', '-V7.3');










