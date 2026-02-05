function h = LL_plot(corr_mat, config)

    % if config.cond == 2
    %     corr_mat = corr_mat_';
    %     y_lab = ['Lag (ms) ', upper(config.obs_subject)];
    %     x_lab = ['Lag (ms) ', upper(config.active_subject)];
    % else
    %     corr_mat = corr_mat_;
    %     y_lab = ['Lag (ms) ', upper(config.active_subject)];
    %     x_lab = ['Lag (ms) ', upper(config.obs_subject)];
    % end
    y_lab = ['Lag (ms) ', upper(config.y_subject)];
    x_lab = ['Lag (ms) ', upper(config.x_subject)];

    n_lags = size(corr_mat,1);

    Lag_X = linspace(0, config.time_length, n_lags);
    Lag_Y = linspace(0, config.time_length, n_lags);

    h.fig = figure;
    h.ax  = axes(h.fig);

    h.im = imagesc(h.ax, Lag_X, Lag_Y, corr_mat);
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