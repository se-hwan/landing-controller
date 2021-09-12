function [X_dn, U_dn, jpos_dn] = data_denormalization_noGRF(nn_data)
    %% parameters
    N = 21; 
    X_tmp = zeros(12, N);
    U_nn_temp = zeros(12, N-1);
    U_tmp = zeros(24, N-1);
    jpos_tmp = zeros(12, N-1);
    dt_val = [0.05 repmat(0.02, 1, 15) [0.05 0.05 0.1 0.5]];
    t_star = zeros(1,N);
    for k = 2:N
        t_star(k) = t_star(k-1) + dt_val(1,k-1);
    end

    load('data/landing_norm_param.mat', 'data_stats');        % data statistics
    
    %% parse data
    X_nn      = reshape(nn_data(1:numel(X_tmp)),size(X_tmp));
    U_nn      = reshape(nn_data(numel(X_tmp) + 1:numel(X_tmp) + numel(U_nn_temp)), size(U_nn_temp));
    jpos_nn   = reshape(nn_data(numel(X_tmp) + numel(U_nn_temp) + 1:numel(X_tmp) + numel(jpos_tmp) + numel(U_nn_temp)), size(jpos_tmp));
    
    %% denormalize data
    X_dn        = X_nn.*data_stats.std_X + data_stats.mean_X;
    U_dn        = U_nn(1:12, :).*data_stats.std_U(1:12, :) + data_stats.mean_U(1:12, :);
    U_dn(13:24, :) = zeros(12, N-1);
    jpos_dn     = jpos_nn.*data_stats.std_jpos + data_stats.mean_jpos;
%     td_dn       = int8(td_nn);
    
    
end
