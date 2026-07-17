
function stats = LL_directional_stats(c_m, mode)

    % Inizializzo sempre
    stats = struct;

    % Costruzione triangoli coerente con LL_plot_detail
    S2K_triangle = triu(c_m, 1);   % S anticipa K
    K2S_triangle = tril(c_m, -1);  % K anticipa S
    diag_values  = diag(c_m);      % sincronia

    eps_val = 1e-6;

    S2K_triangle(abs(S2K_triangle) < eps_val) = NaN;
    K2S_triangle(abs(K2S_triangle) < eps_val) = NaN;
    diag_values(abs(diag_values)   < eps_val) = NaN;

    switch lower(mode)

        case 's2k'
            stats.mean_S2K = mean(S2K_triangle(:), 'omitnan');

        case 'k2s'
            stats.mean_K2S = mean(K2S_triangle(:), 'omitnan');

        case 'diag'
            stats.mean_diag = mean(diag_values, 'omitnan');

        case 'all'
            stats.mean_S2K = mean(S2K_triangle(:), 'omitnan');
            stats.mean_K2S = mean(K2S_triangle(:), 'omitnan');
            stats.mean_diag = mean(diag_values, 'omitnan');

        otherwise
            error('Unknown mode');
    end

end