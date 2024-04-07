function  h=plot_scatterGradient(data_trials,par)
%function h=plot_scatterGradient(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig    =figure;
    par.hfig=hfig;
else
    hfig    =par.hfig;
end 
ms                  = par.ms;
lats                = par.lats;
fine                = par.fine;
widthsub            = par.widthsub;
reverse             = par.reverse;
cmaps               = par.cmaps;
cmapslight          = par.cmapslight;
InField             = par.InField;
InGradient          = par.InGradient;
X_data              = [data_trials.(InField)];  
trialType           = [data_trials.repTrialType];
trialName           = [data_trials.repTrialName];
gradientInfo        = [data_trials.(InGradient)];
nChannels           = length(lats);
% hfig=figure; 
    
% sort input data wrt gradientInfo   
[sortGradients,inds]= sort(gradientInfo);
sortEmb             = X_data(:,inds);
sortName            = trialType(inds);
[classes,idc]       = unique(trialType);
names               = trialName(idc);
nClasses            = length(classes);
subplot(1,widthsub,1:widthsub-nClasses);
hold on; box on; grid on;

uniquePoints        = nan(nClasses,1);
for ic = 1:nClasses
    uniquePoints(ic)=length(unique(sortGradients(sortName==classes(ic))));
end
nColors             = mean(uniquePoints);
ulab                = cell(nClasses,1);
gradientMaps        = cell(nClasses,1);
for ic=1:nClasses    
    sortLlab            = sortName==classes(ic);
    slableft            = sortGradients(sortLlab);
    [~,ulab{ic},lbin]   = histcounts(slableft,nColors);
    colStart            = cmapslight(ic,:);
    colEnd              = cmaps(ic,:);
    gmap                = [linspace(colStart(1),colEnd(1),nColors)' ...
                           linspace(colStart(2),colEnd(2),nColors)' ...
                           linspace(colStart(3),colEnd(3),nColors)'];
    gradientMaps{ic}    = gmap;
    if nChannels==3
        scatter3(sortEmb(lats(1),sortLlab), sortEmb(lats(2),sortLlab), sortEmb(lats(3),sortLlab), ms, gmap(lbin,:), 'filled');
        zlabel(['x' num2str(lats(3))]);
        view(3);
    elseif nChannels==2
        scatter(sortEmb(lats(1),sortLlab), sortEmb(lats(2),sortLlab), ms, gmap(lbin,:), 'filled');
        reverse(3)=false;
    end
    %% sanity check extremalia -> to be deleted
    ifcheck=false;
    if ifcheck
        lwo                 = find(sortLlab);
        [~,wo]              = min(slableft(:,:));
        plot3(sortEmb(lats(1),lwo(wo)), sortEmb(lats(2),lwo(wo)), sortEmb(lats(3),lwo(wo)),'o','markersize',30);
        text (sortEmb(lats(1),lwo(wo)), sortEmb(lats(2),lwo(wo)), sortEmb(lats(3),lwo(wo)),num2str(slableft(wo)),'fontsize',40)
        [~,wo]              = max(slableft(:,:));
        plot3(sortEmb(lats(1),lwo(wo)), sortEmb(lats(2),lwo(wo)), sortEmb(lats(3),lwo(wo)),'o','markersize',30);
        text (sortEmb(lats(1),lwo(wo)), sortEmb(lats(2),lwo(wo)), sortEmb(lats(3),lwo(wo)),num2str(slableft(wo)),'fontsize',40)    
    end
end
xlabel(['x' num2str(lats(1))]);
ylabel(['x' num2str(lats(2))]);
if reverse(1)
    set(gca, 'XDir','reverse')
end
if reverse(2)
    set(gca, 'YDir','reverse');
end
if reverse(3)
    set(gca, 'ZDir','reverse');
end
set(gca,'fontsize',par.fs);


% plot colorbars in different subplots
for ic=1:nClasses
    subplot(1,widthsub,widthsub-nClasses+ic);
    axc = gca;
    colormap(axc,gradientMaps{ic});
    hc  = colorbar;
    % double the weight of color bar
    pos=get(hc,'Position');
    pos(3)=2*pos(3);
    set(hc,'Position',pos)  
    axis off
    ulableft            = ulab{ic};
    uTimes              = size(ulableft,2);
    newTicks            = interp1(linspace(hc.Ticks(1),hc.Ticks(end),uTimes),ulableft, linspace(hc.Ticks(1),hc.Ticks(end),length(hc.Ticks)));
    hc.TickLabels       = round(fine*newTicks)/fine;
    if ic==nClasses
        hc.Label.String     = par.label;
    end
    hc.Title.String     = names{ic};
    set(axc,'fontsize',par.fs);
end   
subplot(1,widthsub,1:widthsub-nClasses);
set(hfig,'Position',par.p);
if nargout>0
    h=hfig;
end
%% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end