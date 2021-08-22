function [X_dn, U_dn, jpos_dn] = data_denormalization(nn_data)
    %% parameters
    N = 21; 
    X_tmp = zeros(12, N);
    U_tmp = zeros(24, N-1);
    jpos_tmp = zeros(12, N-1);
%     td_tmp = zeros(4,1);
    tdx_val = zeros(4, N-1);
    tdx_start = zeros(4, 1);
    dt_val = [0.05 repmat(0.02, 1, 15) [0.05 0.05 0.1 0.5]];
    t_star = zeros(1,N);
    for k = 2:N
        t_star(k) = t_star(k-1) + dt_val(1,k-1);
    end

    load('data/data_stats.mat', 'data_stats');        % data statistics
    
    %% parse data
    X_nn      = reshape(nn_data(1:numel(X_tmp)),size(X_tmp));
    U_nn      = reshape(nn_data(numel(X_tmp) + 1:numel(X_tmp) + numel(U_tmp)), size(U_tmp));
    jpos_nn   = reshape(nn_data(numel(X_tmp) + numel(U_tmp) + 1:numel(X_tmp) + numel(jpos_tmp) + numel(U_tmp)), size(jpos_tmp));
%     td_nn     = reshape(nn_data(numel(X_tmp) + numel(U_tmp) + numel(jpos_tmp) + 1:end), size(td_tmp));
    tdx_val = reshape(nn_data(end-numel(tdx_val)+1:end), N-1, 4).';
    for leg=1:4
        [~, tdx_start(leg)] = max(tdx_val(leg, :));
    end
    f_nn      = U_nn(13:24, :);
    
    %% denormalize data
    X_dn        = X_nn.*data_stats.std_X + data_stats.mean_X;
    U_dn        = U_nn(1:12, :).*data_stats.std_U(1:12, :) + data_stats.mean_U(1:12, :);
    jpos_dn     = jpos_nn.*data_stats.std_jpos + data_stats.mean_jpos;
%     td_dn       = int8(td_nn);

    %% ground reaction force offset
    for leg = 1:4
        xyz_idx = 3*leg-2 : 3*leg;
        f_offset = f_nn(xyz_idx, :);
        f_offset = [zeros(3, tdx_start(leg)-1) f_offset(:, 1:end-tdx_start(leg)+1)];
        
        U_dn(12 + xyz_idx, :) = f_offset.*(data_stats.mass*9.81);
    end
    
end
