function h = LL_plot(corr_mat_, config)

    if config.cond == 2
        corr_mat = corr_mat_';
        y_lab = ['Lag (ms) ', upper(config.obs_subject)];
        x_lab = ['Lag (ms) ', upper(config.active_subject)];
    else
        corr_mat = corr_mat_;
        y_lab = ['Lag (ms) ', upper(config.active_subject)];
        x_lab = ['Lag (ms) ', upper(config.obs_subject)];
    end

    n_lags = size(corr_mat,1);

    Lag_A = linspace(0, 400, n_lags);
    Lag_B = linspace(0, 400, n_lags);

    h.fig = figure;
    h.ax  = axes(h.fig);

    h.im = imagesc(h.ax, Lag_B, Lag_A, corr_mat);
    set(h.ax, 'YDir', 'normal');
    colormap(redbluecmap);
    colorbar;
    caxis([-1 1]);

    xlabel(x_lab);
    ylabel(y_lab);

    title(sprintf('\nCondition %d: \nDirection %d',config.cond, config.dir));

    hold on
    xlims = xlim;
    ylims = ylim;
    plot(xlims, ylims, '--', 'Color',[.5 .5 .5], 'LineWidth',1.2)
    hold off

end