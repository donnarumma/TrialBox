function Y_hat = decoders(lag_X, lag_Y, config)

assert(isfield(config, 'decoder'), 'Missing config.decoder');
assert(isfield(config.decoder, 'method'), 'Missing config.decoder.method');
method = config.decoder.method;


switch lower(method)

    case 'ridge'
        lambdas = config.decoder.ridge_lambdas;
        n_lambda = numel(lambdas);

        n_lags = numel(lag_X);
        [~, n_ch] = size(lag_Y);

        Y_hat = cell(n_lags, n_ch);

        for b = 1:n_lags
            X_ = lag_X{b};
            d  = size(X_,2);
            I  = eye(d);

            % centro X una volta
            Xc = X_ - mean(X_);

            for e = 1:n_ch
                y_ = lag_Y{b,e};
                yc = y_ - mean(y_);

                % GCV 
                gcv_scores = zeros(n_lambda,1);
                for i = 1:n_lambda
                    lambda = lambdas(i);
                    beta   = (Xc' * Xc + lambda * I) \ (Xc' * yc);
                    y_hat_ = Xc * beta;

                    H_diag = sum((Xc / (Xc' * Xc + lambda * I)) .* Xc, 2);
                    gcv_scores(i) = mean(((yc - y_hat_) ./ (1 - H_diag)).^2);
                end

                [~, best_i] = min(gcv_scores);
                lambda_best = lambdas(best_i);
                % Fit finale
                beta = (Xc' * Xc + lambda_best * I) \ (Xc' * yc);
                intercept = mean(y_) - mean(X_) * beta;

                Y_hat{b,e} = X_ * beta + intercept;
            end
        end
   case 'knn'

        assert(isfield(config.decoder, 'knn_k'), ...
            'Missing config.decoder.knn_k');
    
        k = config.decoder.knn_k;
    
        n_lags = numel(lag_X);
        [~, n_ch] = size(lag_Y);
    
        Y_hat = cell(n_lags, n_ch);
    
        for b = 1:n_lags
            X_ = lag_X{b};
    
            %  standardizzazione 
            Xs = zscore(X_);
    
            % precompute neighbors (escludendo self)
            idx = knnsearch(Xs, Xs, 'K', k+1);
            idx = idx(:, 2:end);  % tolgo self
    
            for e = 1:n_ch
                y_ = lag_Y{b,e};
    
                % media sui vicini
                y_hat = mean(y_(idx), 2);
    
                Y_hat{b,e} = y_hat;
            end
        end


    otherwise
        error('Unknown decoder method: %s', method);
end
end
