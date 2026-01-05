function LL_matrix = compute_LL_matrix(A_struct, B_struct, config)

% defaults

assert(config.block_size < config.trial_length, ...
    'block_size must be smaller than trial_length');
% Consistency between corr_obj and use_neural
    switch config.corr_obj
        case "manifold"
            % use_neural can be true or false (future-proof)
        case {"hat-obs", "hat-hat"}
            assert(config.use_neural == true, ...
           'use_neural must be true for corr_obj = %s', config.corr_obj);
        otherwise
            error("Unknown corr_obj");
    end

% filtering data
[A_sel, B_sel]=filter_data(A_struct, B_struct, config)

% build lag structures
if config.use_neural
    [lag_manif_a, lag_manif_b, lag_n_a, lag_n_b] = ...
        convert_to_lag_struct(A_sel, B_sel, config);
else
    [lag_manif_a, lag_manif_b] = ...
        convert_to_lag_struct(A_sel, B_sel, config);
end

%%  select objects to correlatie representation 
switch config.corr_obj

    case "manifold"
        % filter data according to direction and condition
        LL_matrix=f_corr(lag_manif_a, lag_manif_b)
        LL_plot(LL_matrix, config)

    case "hat-obs"
        % Predizioni
        Y_hat_A = decoders(lag_manif_a, lag_n_a, config);
        Y_hat_B = decoders(lag_manif_b, lag_n_b, config);
        
        % Direzione A --> B
        LL_AB = f_corr(Y_hat_A, lag_n_b);
        
        % Direzione B --> A
        LL_BA = f_corr(Y_hat_B, lag_n_a);
        
        % Media simmetrica (come fanno Hasson e Zada nel paper)
        LL_matrix = 0.5 * (LL_AB + LL_BA);
        LL_plot(LL_matrix, config) 

    case "hat-hat" 
       Y_hat_A = decoders(lag_manif_a, lag_n_a, config);
       Y_hat_B = decoders(lag_manif_b, lag_n_b, config); 
       LL_matrix =f_corr(Y_hat_A, Y_hat_B);
       LL_plot(LL_matrix,config)

    otherwise
        error("Unknown corr_obj");
end



end