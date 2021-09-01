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

%% load data
load('data/training_data_landing.mat')

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

%% data normalization

training_data_normalized = struct();

mean_input  = mean(training_data.input, 2);
std_input    = std(training_data.input, 0, 2);

mean_output = mean(training_data.output, 2);
std_output    = std(training_data.output, 0, 2);

mean_X      = reshape(mean_output(1:numel(X_tmp)),size(X_tmp));
mean_U      = reshape(mean_output(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
mean_jpos   = reshape(mean_output(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(jpos_tmp) + numel(U_tmp)), size(jpos_tmp));

std_X       = reshape(std_output(1:numel(X_tmp)),size(X_tmp));
std_U       = reshape(std_output(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
std_jpos    = reshape(std_output(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(jpos_tmp) + numel(U_tmp)), size(jpos_tmp));

data_stats.mean_input = mean_input; data_stats.std_input = std_input;
data_stats.mean_X = mean_X; data_stats.mean_U = mean_U; data_stats.mean_jpos = mean_jpos;
data_stats.std_X = std_X;   data_stats.std_U = std_U;   data_stats.std_jpos = std_jpos;
data_stats.td_scale = 1;
data_stats.mass = mass_val;

save('data/data_stats.mat', 'data_stats')

% prepare normalization values for export

mean_U_offset = mean_U; mean_U_offset(13:24, :) = 0.0;
std_U_offset = std_U; std_U_offset(13:24, :) = mass_val*9.81;
landing_norm_param = [data_stats.mean_input(:); data_stats.std_input(:); mean_X(:); mean_U_offset(:); mean_jpos(:);
                      std_X(:); std_U_offset(:); std_jpos(:)];

fid = fopen('landing_norm_param.dat','w');
fwrite(fid,landing_norm_param,'single');

% normalize input data
training_data_normalized.input = (training_data.input - mean_input)./std_input;
training_data_normalized.output = zeros(length(training_data.output(:,1)) + 4*(N-1), 1);

for data_entry = 1:length(training_data.input)
    %% reshape data
    output_k = training_data.output(:, data_entry);
    
    X_k     = reshape(output_k(1:numel(X_tmp)),size(X_tmp));
    U_k     = reshape(output_k(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
    jpos_k  = reshape(output_k(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(jpos_tmp) + numel(U_tmp)), size(jpos_tmp));

    f_k     = U_k(13:24, :);
%     td      = zeros(4, 1);
    indexed_td = zeros(4, 20);
    
    X_norm = X_tmp; U_norm = U_tmp; jpos_norm = jpos_tmp;
    
    %% ground reaction force normalization
    for leg = 1:4
        % offset ground reaction forces
        xyz_idx = 3*leg-2 : 3*leg;
        f_leg = f_k(xyz_idx, :);
        td_idx = find(f_leg(3,:) > 1, 1);           % how should we normalize this? shouldn't ever be greater than 7-8
        f_leg_offset = f_leg(:, td_idx:end);
        f_leg_offset = [f_leg_offset repmat(f_leg_offset(:, end), 1, td_idx-1)];
        
%         td(leg) = td_idx;                           % store touchdown indices
        indexed_td(leg, td_idx) = 1;
        
        % normalize offset ground reaction forces with bodyweight (could also try min-max normalization)
        min_grf     = min(f_leg_offset,[],2);
        range_grf   = range(f_leg_offset,2);
        f_leg_norm  = (f_leg_offset)./(mass_val*9.81);
        
        U_norm(12 + xyz_idx, :) = f_leg_norm;       % store normalized forces
    end

    %% state normalization
    X_norm = (X_k - mean_X)./std_X; X_norm(1:2, 1) = 0;
    
    %% foot pos. normalization
    U_norm(1:12, :) = (U_k(1:12, :) - mean_U(1:12, :))./std_U(1:12, :);
    
    %% joint pos. normalization
    jpos_norm = (jpos_k - mean_jpos)./std_jpos;

    %% store normalized data
    training_data_normalized.output(:, data_entry) = [X_norm(:); U_norm(:); jpos_norm(:); reshape(indexed_td.', numel(indexed_td), 1)];

end

save('data/training_data_normalized_idx.mat', 'training_data_normalized', '-V7.3')

%% test denormalization
sample_idx = randi(1000, 1);

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

normalized_data_sample = training_data_normalized.output(:, sample_idx);
[X_dn, U_dn, jpos_dn] = data_denormalization(normalized_data_sample);

q_star(1:6,:) = X_dn(1:6,:); qd_star(1:6,:) = X_dn(7:12,:);
q_star(7:18,1:end-1) = jpos_dn; q_star(7:18, end) = q_star(7:18, end-1);

f_star = U_dn(13:24, :); p_star = U_dn(1:12, :);

if show_animation
    showmotion_floatingBase(model,t_star(1:end-1),q_star(:,1:end-1))
end
if show_plots
    plot_results(model, params, t_star, X_dn, U_dn, jpos_dn);
end



%% test sample trajectory
test_idx = sample_idx;

test_input = training_data.input(:, test_idx);
test_trajectory = training_data.output(:, test_idx);

X_star = reshape(test_trajectory(1:numel(X_tmp)),size(X_tmp));
U_star = reshape(test_trajectory(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
jpos_star = reshape(test_trajectory(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(jpos_tmp) + numel(U_tmp)), size(jpos_tmp));

q_star(1:6,:) = X_star(1:6,:); qd_star(1:6,:) = X_star(7:12,:);
q_star(7:18,1:end-1) = jpos_star; q_star(7:18, end) = q_star(7:18, end-1);

f_star = U_star(13:24, :); p_star = U_star(1:12, :);


if show_animation
    showmotion_floatingBase(model,t_star(1:end-1),q_star(:,1:end-1))
end
if show_plots
    plot_results(model, params, t_star, X_star, U_star, jpos_star);
end

