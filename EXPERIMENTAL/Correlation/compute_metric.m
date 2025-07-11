function metric_val = compute_metric(y_test, y_hat, y_train)
    
    y_test = y_test(:);
    y_hat = y_hat(:);
    if nargin < 3
        y_train = [];
    end

    n = length(y_test);

    % Caso generale: almeno due punti e varianze non nulle 
    % okkio che con due punti la correlazione è 1 o -1
    if n > 1 && std(y_test) > 0 && std(y_hat) > 0
        metric_val = corr(y_test, y_hat, 'rows', 'complete');
        return;

    % Caso speciale: più punti ma una varianza nulla (es. y_hat costante) 
    elseif n > 1 && std(y_test) > 0
        SS_res = sum((y_test - y_hat).^2);
        SS_tot = sum((y_test - mean(y_test)).^2);
        r2 = 1 - SS_res / SS_tot;
        r_sign = sign((y_test - mean(y_test))' * (y_hat - mean(y_hat)));
        metric_val = r_sign * r2;
        return;

    % Caso limite: un solo punto (loocv)
   % Caso limite: un solo punto (loocv)
    elseif n == 1
        if ~isempty(y_train) && std(y_train) > 0
            scale = std(y_train);
        elseif std(y_hat) > 0
            scale = std(y_hat);
        else
            % caso estremo max abs val tra y_test y_hat e 1
            scale = max(abs([y_test; y_hat; 1])); % mai 0
        end

    err = abs(y_test - y_hat);
    %% tentativo di imitare la correlazione 
    % errore 0 ottengo 1, se è metà della scala ottengo 0 se è uguale alla 
    % to dà 1 se previsione perfetta
    % se err = scale/2 ottengo errore medio metric =0
    % se err = scale erroe alto metric=-1
    % da invertire
    metric_val = 1 - 2 * err / scale;
    % clamp tra -1 e 1
    metric_val = max(-1, min(1, metric_val)); 
    metric_val = double(metric_val); % <-- forza scalare
    return;

    % Caso degenere
    else
        metric_val = 0;
        return;
    end
end
