function   h = plot_scatter(data_trials,par)
% function h = plot_scatter(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig    =figure;
    par.hfig=hfig;
else
    hfig    =par.hfig;
end 
hold on; box on; grid on;
set(hfig,'defaultLegendAutoUpdate','off');
xfld            = par.xfld;
evaltime        = par.evaltime;
InField         = par.InField;
nTrials         = length(data_trials);
markerSize      = par.markerSize;
marker          = par.marker;
nChannels       = size(data_trials(1).(InField),1);
if isempty(par.wd)
    idxComponents= 1:nChannels;
else
    idxComponents= par.wd;      
end
if length(idxComponents)>3
    idxComponents=idxComponents(1:3);
end
if length(idxComponents)<2
    return
end

X_data          =  nan(length(data_trials),length(idxComponents));
for iTrial=1:nTrials
    timeField   = data_trials(iTrial).([xfld InField]);
    if isempty(evaltime)
        iTime=1;
    else
        [~,iTime]   = min(abs(timeField-evaltime));
    end
    X_data(iTrial,:) = data_trials(iTrial).(InField)(idxComponents,iTime)';
end
trialTypes      = [data_trials.trialType];
trialNames      = {data_trials.trialName};
[types,idx]     = unique(trialTypes);
labels          = trialNames(idx);
nClasses        = length(types);
InField_title   = strcat(InField, ' scatter plot');
% markerSize      = 15;
if isempty(par.cmaps)
    colors      = jet(nClasses);
else
    colors      = par.cmaps;
end
tickLabels      = labels; 
nComponents     = length(idxComponents);
tstart          = 1/(2*nClasses);
ticks           = linspace(0,1,nClasses+1)+tstart;
ticks           = ticks(1:end-1);
filled          = par.filled;
notLegend       = true;
% 'x','markersize',ms,'MarkerEdgeColor',col,
hfigs           = gobjects(nTrials,1); 
% if nComponents == 2
    % scatter (X_data(:,1),X_data(:,2),            markerSize,colors(trialTypes,:),'filled')
% elseif nComponents==3
    % scatter3(X_data(:,1),X_data(:,2),X_data(:,3),markerSize,colors(trialTypes,:),'filled','x')
    % plot3   (X_data(:,1),X_data(:,2), X_data(:,3),'LineStyle','none','marker','x','markersize',markerSize,'MarkerEdgeColor',colors(trialTypes,:),'MarkerFaceColor',colors(trialTypes,:),'linewidth',  5);
    colfaces = colors(trialTypes,:);
    if filled
        coledges=colfaces;
    else 
        coledges=repmat({'k'},nTrials,1);
    end
        
    for iTrial=1:nTrials
        if nComponents==2
            hfigs(iTrial)=plot   (X_data(iTrial,1),X_data(iTrial,2),    'LineStyle','none', ...
                                                                        'marker',marker,'MarkerSize',markerSize,...
                                                                        'MarkerEdgeColor',coledges(iTrial,:), ...
                                                                        'MarkerFaceColor',colfaces(iTrial,:));

        else
            hfigs(iTrial)=plot3   (X_data(iTrial,1),X_data(iTrial,2), X_data(iTrial,3),'LineStyle','none', ...
                                                                        'marker',marker,'MarkerSize',markerSize,...
                                                                        'MarkerEdgeColor',coledges(iTrial,:), ...
                                                                        'MarkerFaceColor',colfaces(iTrial,:));
        end
        if notLegend
            hfigs(iTrial).Annotation.LegendInformation.IconDisplayStyle='off';
        end
                                                                    % 'linewidth',  2.5);
    end
    % scatter3   (X_data(:,1),X_data(:,2), X_data(:,3), markerSize, repmat({'MarkerFaceColor'},nTrials,1),colors(trialTypes,:),repmat({'MarkerEdgeColor'},nTrials,1),colors(trialTypes,:),marker); 
% end

title(InField_title)
colormap(colors);
c              = colorbar;
c.Label.String = 'Labels';
set(c, 'Ticks', ticks, 'TickLabels', tickLabels);

if nargout>0
    h=hfig;
end
%% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end