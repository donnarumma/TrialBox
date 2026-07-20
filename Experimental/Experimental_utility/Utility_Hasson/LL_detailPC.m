function [c_m, stats] = LL_detailPC(corr_mat, config, t_start_y,t_end_y, t_start_x, t_end_x)
    
    % ------------------------------------------------------------
    % INPUT:
    % corr_mat : matrice lag-lag completa [n_lags x n_lags]
    % quadrata e uniformemente campionata nel tempo
    % config   : struttura con t_start, t_end, soggetti ecc.
    % eventuali limiti temporali per ritagliare la matrice
    %
    % OUTPUT:
    % c_m   : sottomatrice ritagliata
    % h     : handle grafico (figura, assi, image)
    % stats : (opzionale) statistiche direzionali
    % ------------------------------------------------------------

    n_lags = size(corr_mat,1);
    full_length = config.t_end - config.t_start;

    t_start_x = max(config.t_start, min(config.t_end, t_start_x));
    t_end_x   = max(config.t_start, min(config.t_end, t_end_x));
    t_start_y = max(config.t_start, min(config.t_end, t_start_y));
    t_end_y   = max(config.t_start, min(config.t_end, t_end_y));

    % map time to matrix index (n.b.! assume a linear relationship between
    % mapping and time)
    t_start_x_ = floor(1 + (t_start_x - config.t_start) / full_length * (n_lags - 1));
    t_end_x_   = ceil(1 + (t_end_x  - config.t_start) / full_length * (n_lags - 1));
    t_start_x_= max(1, min(n_lags, t_start_x_));
    t_end_x_  = max(1, min(n_lags, t_end_x_));
    t_start_y_ = floor(1 + (t_start_y - config.t_start) / full_length * (n_lags - 1));
    t_end_y_   = ceil(1 + (t_end_y  - config.t_start) / full_length * (n_lags - 1));
    t_start_y_= max(1, min(n_lags, t_start_y_));
    t_end_y_  = max(1, min(n_lags, t_end_y_));
   

    idx_x= t_start_x_:t_end_x_;
    idx_y= t_start_y_:t_end_y_;
 
    c_m=corr_mat(idx_y, idx_x);

    % ------------------------------------------------------------
    % Calcolo statistiche solo se richiesto
    % ------------------------------------------------------------
    stats = LL_directional_stats(c_m, 'all');
end