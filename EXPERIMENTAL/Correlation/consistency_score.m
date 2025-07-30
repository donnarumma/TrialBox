function R2 = consistency_score(X, Y)
    % p is likely equal to m
    % n × p
    X = zscore(X);  
     % n × m
    Y = zscore(Y); 

    % Fit multivariate linear model: solve B in X*B ≈ Y
    B = X \ Y;          
    Y_hat = X * B;

    % Compute residual sum of squares (SSE)
    residual = Y - Y_hat;
    SSE = sum(residual.^2, 'all');

    % Compute total sum of squares (SST)
    SST = sum((Y - mean(Y)).^2, 'all');

    % Compute R^2
    R2 = 1 - (SSE / SST);
end
