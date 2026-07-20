function plotHassonProfiles( ...
    S2K_time,...
    K2S_time,...
    t_vec,...
    condName,...
    RT_mean,...
    ExT_mean,...
    save_dir,...
    compName)
% PLOTHASSONPROFILES
%
% Purpose
%   Generates temporal visualizations of Hasson directional coupling
%   profiles between monkey S and monkey K neural manifolds.
%
%   The function produces single-direction, multi-direction and
%   condition-averaged summary plots for both S->K and K->S coupling
%   trajectories.
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
%   PNG, PDF and MAT figures saved to disk.
%
% Used by
%   MAIN_HASSON_WITH_EXTRACT_SPIKES
%
% Author
%   Mirco Frosolone
%

%% Limits

maxValue = max([S2K_time(:);K2S_time(:)]) * 1.1;
minValue = min([S2K_time(:);K2S_time(:)]) * 1.1;

%% Output directory

save_dirLINplot = fullfile(save_dir,'LINEAR_PLOT',compName);

%% ============================================================
%  Single direction plots
%  ============================================================

for ncond = 1:3

    for ndir = 1:8

        figure

        plot(t_vec,squeeze(S2K_time(ncond,ndir,:)), ...
            'b','LineWidth',2)

        hold on

        plot(t_vec,squeeze(K2S_time(ncond,ndir,:)), ...
            'g','LineWidth',2)

        xline(RT_mean,'r--','RT','LineWidth',2,...
            'LabelHorizontalAlignment','center');

        xline(ExT_mean,'m:','ExT','LineWidth',2,...
            'LabelHorizontalAlignment','center');

        yline(0,'k-','LineWidth',1,...
            'HandleVisibility','off');

        xlabel('Time (ms)')
        ylabel('Hasson value')
        ylim([minValue maxValue])

        legend('S2K','K2S')

        title([condName{ncond} ...
            '\newline      Dir ' num2str(ndir)])

        grid on

        set(gcf,'WindowState','maximized');
        set(gcf,'Units','normalized',...
            'OuterPosition',[0 0 1 1]);

        par.savePlotEpsPdfMat.dir_png = ...
            fullfile(save_dirLINplot,'PNGs');

        par.savePlotEpsPdfMat.dir_pdf = ...
            fullfile(save_dirLINplot,'PDFs');

        par.savePlotEpsPdfMat.dir_mat = ...
            fullfile(save_dirLINplot,'MATfiles');

        par.savePlotEpsPdfMat.file_name = ...
            ['Sing_',condName{ncond},'_Dir',num2str(ndir)];

        savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)

        close(gcf)

    end

end

%% ============================================================
%  All directions subplot
%  ============================================================

for ncond = 1:3

    figure

    for ndir = 1:8

        subplot(2,4,ndir)

        plot(t_vec,squeeze(S2K_time(ncond,ndir,:)), ...
            'b','LineWidth',2)

        hold on

        plot(t_vec,squeeze(K2S_time(ncond,ndir,:)), ...
            'g','LineWidth',2)

        xline(RT_mean,'r--','RT','LineWidth',2,...
            'LabelHorizontalAlignment','center');

        xline(ExT_mean,'m:','ExT','LineWidth',2,...
            'LabelHorizontalAlignment','center');

        yline(0,'k-','LineWidth',1,...
            'HandleVisibility','off');

        title(['Dir ' num2str(ndir)])

        xlabel('Time [ms]')
        ylabel('Hasson value')

        ylim([minValue maxValue])

        grid on

    end

    axes('Position',[0 0 1 1],'Visible','off');

    line(NaN,NaN,'LineWidth',2,'Color','b');
    line(NaN,NaN,'LineWidth',2,'Color','g');

    legend({'S2K','K2S'},...
        'Location','south');

    sgtitle(condName{ncond})

    set(gcf,'WindowState','maximized');
    set(gcf,'Units','normalized',...
        'OuterPosition',[0 0 1 1]);

    par.savePlotEpsPdfMat.dir_png = ...
        fullfile(save_dirLINplot,'PNGs');

    par.savePlotEpsPdfMat.dir_pdf = ...
        fullfile(save_dirLINplot,'PDFs');

    par.savePlotEpsPdfMat.dir_mat = ...
        fullfile(save_dirLINplot,'MATfiles');

    par.savePlotEpsPdfMat.file_name = ...
        ['Sub_',condName{ncond},'_allDir'];

    savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)

    close(gcf)

end

%% ============================================================
%  Mean direction plot
%  ============================================================

for ncond = 1:3

    S2K_mean = squeeze(mean(S2K_time(ncond,:,:),2));
    K2S_mean = squeeze(mean(K2S_time(ncond,:,:),2));

    figure

    plot(t_vec,S2K_mean,...
        'b','LineWidth',3)

    hold on

    plot(t_vec,K2S_mean,...
        'g','LineWidth',3)

    xline(RT_mean,'r--','RT','LineWidth',2,...
        'LabelHorizontalAlignment','center');

    xline(ExT_mean,'m:','ExT','LineWidth',2,...
        'LabelHorizontalAlignment','center');

    yline(0,'k-','LineWidth',1,...
        'HandleVisibility','off');

    xlabel('Time [ms]')
    ylabel('Mean Hasson value')

    ylim([minValue maxValue])

    legend('S2K','K2S')

    title([condName{ncond} ...
        '\newline   ALL Dir'])

    grid on

    set(gcf,'WindowState','maximized');
    set(gcf,'Units','normalized',...
        'OuterPosition',[0 0 1 1]);

    par.savePlotEpsPdfMat.file_name = ...
        ['MEAN_',condName{ncond},'_allDir'];

    savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)

    close(gcf)

end

%% ============================================================
%  Mean direction subplots
%  ============================================================

figure

for ncond = 1:3

    S2K_mean = squeeze(mean(S2K_time(ncond,:,:),2));
    K2S_mean = squeeze(mean(K2S_time(ncond,:,:),2));

    subplot(3,1,ncond)

    plot(t_vec,S2K_mean,...
        'b','LineWidth',3)

    hold on

    plot(t_vec,K2S_mean,...
        'g','LineWidth',3)

    xline(RT_mean,'r--','RT','LineWidth',2,...
        'LabelHorizontalAlignment','center');

    xline(ExT_mean,'m:','ExT','LineWidth',2,...
        'LabelHorizontalAlignment','center');

    yline(0,'k-','LineWidth',1,...
        'HandleVisibility','off');

    xlabel('Time [ms]')
    ylabel('Mean Hasson value')

    legend('S2K','K2S')

    title([condName{ncond} ...
        '\newline   ALL Dir'])

end

axes('Position',[0 0 1 1],'Visible','off');

line(NaN,NaN,'LineWidth',2,'Color','b');
line(NaN,NaN,'LineWidth',2,'Color','g');

legend({'S2K','K2S'},...
    'Location','south');

set(gcf,'WindowState','maximized');
set(gcf,'Units','normalized',...
    'OuterPosition',[0 0 1 1]);

par.savePlotEpsPdfMat.file_name = ...
    'Sub_MEAN_allDir';

savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)

close(gcf)

end