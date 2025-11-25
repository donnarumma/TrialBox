% function to compute correlation and control limit cases
function metric_val = compute_metric(y_test, y_hat, y_train)
    % y_train rimane per compatibilità, qui non usato.
    if nargin < 3
        y_train = [];
    end

    % Vettori colonna
    y_test = y_test(:);
    y_hat  = y_hat(:);

    % 1) Usa solo coppie valide (pairwise complete).
    %    Se dopo il filtro non resta nulla (nessuna coppa valida)
    % → NaN (non 0).
    valid  = ~isnan(y_test) & ~isnan(y_hat);
    y_test = y_test(valid);
    y_hat  = y_hat(valid);

    n   = length(y_test);
    tol = 1e-12;

    %n <= 2 (r ±1 degenere con due punti)  → NaN
    %    - var(y_test)==0 o var(y_hat)==0 → NaN
    %    - altrimenti: Pearson r 
    if n == 0
        metric_val = NaN;

    elseif n <= 2
        metric_val = NaN;

    else
        sy = std(y_test);
        sh = std(y_hat);

        if sy <= tol || sh <= tol
            metric_val = NaN;

        else
        % Pearson corr (già filtrati i NaN sopra)
            r = corr(y_test, y_hat);   
        
            % Clamp di sicurreza tra -1 e 1. (overshooting possibile)
            if abs(r) > 1 && abs(abs(r) - 1) <= tol
                r = sign(r);
            end
        
            metric_val = r;
        end

    end

    metric_val = double(metric_val);
end
