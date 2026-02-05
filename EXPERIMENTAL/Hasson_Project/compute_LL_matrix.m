function LL_matrix = compute_LL_matrix(X_struct, Y_struct, config)

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
[X_sel, Y_sel]=filter_data(X_struct, Y_struct, config)

% build lag structures
if config.use_neural
    [lag_manif_X, lag_manif_Y, lag_n_X, lag_n_Y] = ...
        convert_to_lag_struct(X_sel, Y_sel, config);
else
    [lag_manif_X, lag_manif_Y] = ...
        convert_to_lag_struct(X_sel, Y_sel, config);
end

%%  select objects to correlatie representation 
switch config.corr_obj

    case "manifold"
        % filter data according to direction and condition
        LL_matrix=f_corr(lag_manif_X, lag_manif_Y)

    case "hat-obs"
        % Predizioni
        Y_hat_X = decoders(lag_manif_X, lag_n_X, config);
        Y_hat_Y = decoders(lag_manif_Y, lag_n_Y, config);
        
        % Direzione A --> B
        LL_XY = f_corr(Y_hat_X, lag_n_Y);
        
        % Direzione B --> A
        LL_YX = f_corr(Y_hat_Y, lag_n_X);
        
        % Media simmetrica (come fanno Hasson e Zada nel paper)
        LL_matrix = 0.5 * (LL_XY + LL_YX);
         

    case "hat-hat" 
       Y_hat_X = decoders(lag_manif_X, lag_n_X, config);
       Y_hat_Y = decoders(lag_manif_Y, lag_n_Y, config); 
       LL_matrix =f_corr(Y_hat_X, Y_hat_Y);

    otherwise
        error("Unknown corr_obj");
end



end