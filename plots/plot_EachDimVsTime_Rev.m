function h = plot_EachDimVsTime_Rev(data_trials, par)
% function h = plot_EachDimVsTime_Rev(data_trials, par)
execinfo = par.exec;
if ~isempty(execinfo); t = tic; fprintf('Function: %s ', mfilename); end

if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig = figure;
    par.hfig = hfig;
else
    hfig = par.hfig;
end

nTrials     = length(data_trials);
keep        = par.keep;
cmaps       = par.cmaps;
cmapslight  = par.cmapslight;
grayzone    = par.grayzone;
nCols       = par.nCols;
InField     = par.InField;

nComponentsAll = size(data_trials(1).(InField),1);
if isfield(par,'keepComponents') && ~isempty(par.keepComponents)
    keepComp = par.keepComponents;
else
    keepComp = 1:nComponentsAll;
end
nComponents = numel(keepComp);

pos = get(hfig, 'Position');
set(hfig, 'Position', [pos(1) pos(2) 2*pos(3) pos(4)]);

xfld        = 'time';
timextkl    = data_trials(1).([xfld InField]);
nTimes      = size(timextkl,2);
nRows       = ceil(nComponents / nCols);
[~, idx]    = unique([data_trials.trialType]);
trialNames  = {data_trials.trialName};
names       = trialNames(idx);
newk        = 1:nComponents;
nRows_mod   = nRows;
nCols_mod   = nCols;

toleg = {};
cc    = gobjects(0);
lw    = 2.5;
ylcs  = cell(nComponents,1);

for iTrial = 1:nTrials
    col  = cmaps(iTrial,:);
    col2 = cmapslight(iTrial,:);
    for iComponent = 1:nComponents
        compIdx = keepComp(iComponent);

        if iTrial == 1
            ylcs{iComponent} = [+inf, -inf];
        end

        subplot(nRows_mod, nCols_mod, newk(iComponent));
        hold on;

        if ismember(iTrial, keep)
            mtl = data_trials(iTrial).(InField)(compIdx,:);
            stl = data_trials(iTrial).([InField 'sigma'])(compIdx,:);

            if par.novariance
                stl = stl*0;
            end

            ylcs{iComponent} = [ ...
                min([min(mtl-stl), ylcs{iComponent}(1)]), ...
                max([max(mtl+stl), ylcs{iComponent}(2)])];

            par.times = timextkl;
            par.c1    = col;
            par.c2    = col2;
            par.c3    = col;
            par.fill  = 1;
            par.hfig  = hfig;

            forP = plotMeanAndVariancePlus(mtl, stl, par);
            plot(timextkl, mtl, 'LineWidth', lw, 'Color', col);

            if iTrial == keep(end)
                try
                    ylmh = par.YLIM;
                catch
                    ylmh = ylcs{iComponent};
                end

                hh   = diff(ylmh);
                dymh = hh/100;

                if length(grayzone) == 2
                    mg = grayzone(1);
                    sg = grayzone(2);
                else
                    mg = mean(grayzone);
                    sg = std(grayzone);
                end

                x = mg - sg;
                y = ylmh(1);
                w = sg * 2;

                try
                    mh1 = rectangle('Position',[x,y,w,hh], ...
                                    'FaceColor',0.9*ones(3,1), ...
                                    'LineStyle','none');
                    uistack(mh1,'bottom');
                catch
                end

                mh2 = plot([mg,mg], ylmh, 'k--');
                uistack(mh2,'bottom');

                if ~par.novariance
                    ylim([ylmh(1)-dymh, ylmh(2)+dymh]);
                end
            end

            if iComponent == 1
                toleg{end+1} = names{iTrial};
                if nTimes > 1
                    cc(end+1) = forP.hmean;
                end
            end
        end
    end
end

selY = newk(1):nCols_mod:max(newk);
selX = newk(max(1,end-nCols+1):end);

for iComponent = 1:nComponents
    compIdx = keepComp(iComponent);

    subplot(nRows_mod, nCols_mod, newk(iComponent));
    box on;
    grid on;

    dx = (timextkl(2)-timextkl(1))/2;
    xl = [timextkl(1)-dx, timextkl(end)+dx];
    xlim(xl);

    try
        title(par.titles{compIdx}, 'Interpreter', 'latex', 'FontSize', 16);
    catch
    end

    if ismember(newk(iComponent), selX) && nTimes > 1
        xlabel('Time (ms)');
    end

    if ismember(newk(iComponent), selY)
        ylabel(par.ylabel, 'Interpreter', 'latex');
    end
end

if isfield(par,'legendLocation') && ~isempty(cc)
    if strcmpi(par.legendLocation,'right')
        ax = findobj(hfig,'Type','Axes');
        ax = ax(~strcmp(get(ax,'Tag'),'legend'));

        for iax = 1:numel(ax)
            posAx = get(ax(iax),'Position');
            posAx(3) = posAx(3) * 0.82;
            set(ax(iax),'Position',posAx);
        end

        hLegend = legend(cc, toleg);
        set(hLegend, 'Units', 'normalized');
        set(hLegend, 'Orientation', 'vertical');
        set(hLegend, 'Box', 'on');
        set(hLegend, 'Color', 'w');
        set(hLegend, 'EdgeColor', [0.8 0.8 0.8]);
        set(hLegend, 'FontSize', 10);
        hLegend.ItemTokenSize = [10, 10];
        set(hLegend, 'Position', [0.85 0.43 0.11 0.20]);

    elseif strcmpi(par.legendLocation,'bottom')
        ax = findobj(hfig,'Type','Axes');
        ax = ax(~strcmp(get(ax,'Tag'),'legend'));

        for iax = 1:numel(ax)
            posAx = get(ax(iax),'Position');
            posAx(2) = posAx(2) + 0.05;
            posAx(4) = posAx(4) * 0.90;
            set(ax(iax),'Position',posAx);
        end

        hLegend = legend(cc, toleg);
        set(hLegend, 'Units', 'normalized');
        set(hLegend, 'Position', [0.25 0.01 0.50 0.05]);
        set(hLegend, 'Orientation', 'horizontal');
        set(hLegend, 'FontSize', 10);
        set(hLegend, 'Box', 'on');

    else
        hLegend = legend(cc, toleg);
        set(hLegend, 'Location', 'best');
    end
end

if nargout > 0
    h = hfig;
end

if ~isempty(execinfo)
    out.exectime = toc(t);
    fprintf('| Time Elapsed: %.2f s\n', out.exectime);
end