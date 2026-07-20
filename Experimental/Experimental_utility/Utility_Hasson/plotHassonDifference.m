function plotHassonDifference( ...
    S2K_time,...
    K2S_time,...
    t_vec,...
    condName,...
    RT_mean,...
    ExT_mean,...
    save_dir,...
    compName)
% PLOTHASSONDIFFERENCE
%
% Purpose
%   Computes and visualizes directional asymmetries in Hasson coupling
%   between monkey S and monkey K.
%
%   The directional difference is defined as:
%
%       K2S - S2K
%
%   Positive values indicate stronger K->S coupling, whereas negative
%   values indicate stronger S->K coupling.
%
% Inputs
%   S2K_time
%       Hasson coupling values from monkey S to monkey K
%       [condition x direction x time]
%
%   K2S_time
%       Hasson coupling values from monkey K to monkey S
%       [condition x direction x time]
%
%   t_vec
%       Time vector used for temporal visualization [ms]
%
%   condName
%       Cell array containing condition labels
%
%   RT_mean
%       Mean behavioural reaction time [ms]
%
%   ExT_mean
%       Mean behavioural exit time [ms]
%
%   save_dir
%       Root directory for figure storage
%
%   compName
%       PCA component-group identifier
%
% Outputs
%   PNG, PDF and MAT figures highlighting K->S versus S->K dominance.
%
% Used by
%   MAIN_HASSON_WITH_EXTRACT_SPIKES
%
% Author
%   Mirco Frosolone
%


%% Difference matrix

Diff_time = K2S_time - S2K_time;

colK2S    = [0.10 0.55 0.25];
colS2K    = [0.15 0.35 0.85];
baseColor = [0.25 0.25 0.25];

minDiffValue = min(Diff_time(:),[],'omitnan');
maxDiffValue = max(Diff_time(:),[],'omitnan');

save_dirLINplot = fullfile( ...
    save_dir,...
    'LINEAR_PLOT_DIFF',...
    compName);

savePlotCfg.dir_png = fullfile(save_dirLINplot,'PNGs');
savePlotCfg.dir_pdf = fullfile(save_dirLINplot,'PDFs');
savePlotCfg.dir_mat = fullfile(save_dirLINplot,'MATfiles');

%% ============================================================
%  Single directions
%  ============================================================

for ncond = 1:3

    for ndir = 1:8

        y = squeeze(Diff_time(ncond,ndir,:))';

        yPos = y;
        yPos(yPos <= 0) = NaN;

        yNeg = y;
        yNeg(yNeg >= 0) = NaN;

        figure
        hold on

        hBase = plot(t_vec,y,...
            '-','Color',baseColor,...
            'LineWidth',1.2);

        hPos = plot(t_vec,yPos,...
            '-','Color',colK2S,...
            'LineWidth',3);

        hNeg = plot(t_vec,yNeg,...
            '-','Color',colS2K,...
            'LineWidth',3);

        xline(RT_mean,'r--','RT','LineWidth',2,...
            'LabelHorizontalAlignment','center');

        xline(ExT_mean,'m:','ExT','LineWidth',2,...
            'LabelHorizontalAlignment','center');

        yline(0,'k-','LineWidth',1,...
            'HandleVisibility','off');

        xlabel('Time [ms]')
        ylabel('Diff Hasson value [K2S - S2K]')

        ylim([minDiffValue maxDiffValue])

        legend([hPos hNeg hBase], ...
            {'K2S > S2K',...
             'S2K > K2S',...
             'K2S - S2K profile'},...
             'Location','best')

        title([condName{ncond} ...
            '\newline      Dir ' num2str(ndir)])

        grid on

        set(gcf,'WindowState','maximized');
        set(gcf,'Units','normalized',...
            'OuterPosition',[0 0 1 1]);

        savePlotCfg.file_name = ...
            ['Sing_DIFF_', ...
            condName{ncond}, ...
            '_Dir', ...
            num2str(ndir)];

        savePlotEpsPdfMat(gcf,savePlotCfg)

        close(gcf)

    end
end

%% ============================================================
%  All directions subplot
%  ============================================================

for ncond = 1:3

    figure

    for ndir = 1:8

        y = squeeze(Diff_time(ncond,ndir,:))';

        yPos = y;
        yPos(yPos <= 0) = NaN;

        yNeg = y;
        yNeg(yNeg >= 0) = NaN;

        subplot(2,4,ndir)
        hold on

        plot(t_vec,y,...
            '-','Color',baseColor,...
            'LineWidth',1.2,...
            'HandleVisibility','off');

        plot(t_vec,yPos,...
            '-','Color',colK2S,...
            'LineWidth',3,...
            'HandleVisibility','off');

        plot(t_vec,yNeg,...
            '-','Color',colS2K,...
            'LineWidth',3,...
            'HandleVisibility','off');

        xline(RT_mean,'r--','RT','LineWidth',2,...
            'LabelHorizontalAlignment','center');

        xline(ExT_mean,'m:','ExT','LineWidth',2,...
            'LabelHorizontalAlignment','center');

        yline(0,'k-','LineWidth',1,...
            'HandleVisibility','off');

        title(['Dir ' num2str(ndir)])

        xlabel('Time [ms]')
        ylabel('Diff Hasson')

        ylim([minDiffValue maxDiffValue])

        grid on

    end

    axes('Position',[0 0 1 1],'Visible','off');

    h1 = line(NaN,NaN,'LineWidth',3,'Color',colK2S);
    h2 = line(NaN,NaN,'LineWidth',3,'Color',colS2K);
    h3 = line(NaN,NaN,'LineWidth',1.2,'Color',baseColor);

    legend([h1 h2 h3], ...
        {'K2S > S2K',...
         'S2K > K2S',...
         'K2S - S2K profile'},...
         'Location','south');

    sgtitle(condName{ncond})

    set(gcf,'WindowState','maximized');
    set(gcf,'Units','normalized',...
        'OuterPosition',[0 0 1 1]);

    savePlotCfg.file_name = ...
        ['Sub_DIFF_', ...
        condName{ncond}, ...
        '_allDir'];

    savePlotEpsPdfMat(gcf,savePlotCfg)

    close(gcf)

end

%% ============================================================
%  Mean direction plot
%  ============================================================

for ncond = 1:3

    Diff_mean = squeeze( ...
        mean(Diff_time(ncond,:,:),2,'omitnan'))';

    yPos = Diff_mean;
    yPos(yPos <= 0) = NaN;

    yNeg = Diff_mean;
    yNeg(yNeg >= 0) = NaN;

    figure
    hold on

    hBase = plot(t_vec,Diff_mean,...
        '-','Color',baseColor,...
        'LineWidth',1.2);

    hPos = plot(t_vec,yPos,...
        '-','Color',colK2S,...
        'LineWidth',3);

    hNeg = plot(t_vec,yNeg,...
        '-','Color',colS2K,...
        'LineWidth',3);

    xline(RT_mean,'r--','RT','LineWidth',2,...
        'LabelHorizontalAlignment','center');

    xline(ExT_mean,'m:','ExT','LineWidth',2,...
        'LabelHorizontalAlignment','center');

    yline(0,'k-','LineWidth',1,...
        'HandleVisibility','off');

    xlabel('Time [ms]')
    ylabel('Mean Diff Hasson value [K2S - S2K]')

    ylim([minDiffValue maxDiffValue])

    legend([hPos hNeg hBase], ...
        {'K2S > S2K',...
         'S2K > K2S',...
         'Mean K2S - S2K'},...
         'Location','best')

    title([condName{ncond} ...
        '\newline   ALL Dir'])

    grid on

    set(gcf,'WindowState','maximized');
    set(gcf,'Units','normalized',...
        'OuterPosition',[0 0 1 1]);

    savePlotCfg.file_name = ...
        ['MEAN_DIFF_', ...
        condName{ncond}, ...
        '_allDir'];

    savePlotEpsPdfMat(gcf,savePlotCfg)

    close(gcf)

end

end