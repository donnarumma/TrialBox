function   h = plot_pValues(p_classes,par)
% function h = plot_pValues(p_classes,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig    =figure;
    par.hfig=hfig;
else
    hfig    =par.hfig;
end
InField         = par.InField;
xfld            = par.xfld;
nCols           = par.nCols;
% nRows           = par.nRows;
time            = p_classes(1).([xfld InField]);
dt              = par.dt;
% dt              = 50; % ms
inds            = time>=time(1)+dt & time<=time(end)-dt;
nClasses        = length(p_classes);
nChannels       = size(p_classes(1).(InField),1);
nRows           = ceil(nChannels / nCols);
explained       = par.explained;
newk            = 1:nChannels;
legplot         = 2;
lw              = 4;
cc              = [];    
if legplot==2
    hSub        = subplot(nRows+1, 1, 1);
    set(gca,'xticklabel','');set(gca,'yticklabel','');
    axis off;
    newk        = indexesForSubplotsWithLegend2(nRows,nCols);
    nRows_mod   = nRows+1;
    nCols_mod   = nCols;
elseif legplot==1
    hSub    = subplot(1, nCols+1, 1);
    set(gca,'xticklabel','');set(gca,'yticklabel','');
    axis off;
    newk=indexesForSubplotsWithLegend(nRows,nCols);
    nRows_mod=nRows;
    nCols_mod=nCols+1;
else
    nRows_mod=nRows;
    nCols_mod=nCols;
end
selY=newk(1):nCols_mod:max(newk);
selX=newk(end-nCols+1):newk(end);
if nClasses>1
    nComps  = nchoosek(nClasses,2);
else
    nComps = 1;
end
% cmaps   = linspecer(nComps);
if nClasses==1
    cmaps   = 0.2*ones(1,3);
else
    cmaps   = linspecer(nClasses+nComps+10);

    cmaps   = cmaps(1:2:end,:);
end
limval  = [inf,-inf];
totit   = cell(1,nChannels);
for ichannel=1:nChannels
    totit{ichannel} = sprintf('$$\\tilde{\\mathbf x}_{%d,t}$$',ichannel);
    if ~isempty(explained)
        totit{ichannel} = sprintf('%s (%2.1f)',totit{ichannel},explained(ichannel));
    end
end
for iChannel=1:nChannels
    subplot(nRows_mod,nCols_mod,newk(iChannel))
    hold on; box on; grid on;
    set(gca,'fontsize',13)
    cc=[];
    toleg=cell(0,0);
    indcol  = 0;
    for idp=1:length(par.decisions)
        % xline(0,'label','Decision Point','LabelVerticalAlignment','bottom','LineWidth',5);
        xline(par.decisions(idp),'label',par.decisionsN{idp},'LabelVerticalAlignment','bottom','LineWidth',5);
    end
    yline(0.05,'LineStyle','--');
    yline(0.1,'LineStyle','--');
    yline(0.5,'LineStyle','--');
    if false
        cc(end+1)=plot(time,pVals(iChannel,:),'Color',cmaps(end,:),'LineWidth',lw);
        trialNames = {p_classes.trialName};
        toleg{end+1}=trialNames{1};
        for itn=2:length(trialNames)
            toleg{end}=sprintf('%s vs %s',toleg{end},trialNames{itn});
        end
        % toleg{end+1}=sprintf('%s vs %s',pComps(iClass).trialName,pComps(ic).trialName);
        limval(1)=min([limval(1);pVals(:)]);
        limval(2)=max([limval(2);pVals(:)]);
    else
        for iClass=1:nClasses
            for ic=1:nClasses
                if ic<=iClass & nClasses > 1
                    continue
                else
                    indcol=indcol+1;
                    % col = cmaps(iClass,:)+cmaps(ic,:);
                    % if any(col>1)
                    %     col=col/max(col);
                    % end
                    % col=col*0.8;
                    % cc(end+1)=plot(time,comp,'Color',col,'LineWidth',lw);
                    comp        = squeeze(p_classes(iClass).(InField)(iChannel,:,ic));
                    comp        = smooth(comp);
                    plfunc      = str2func('plot'); fcolor = 'Color'; xval=time(inds);
                    if length(inds)==1
                        plfunc=str2func('bar');     fcolor='FaceColor'; xval=indcol;
                    end
                    cc(end+1)   = plfunc(xval,comp(inds),fcolor,cmaps(indcol,:),'LineWidth',lw);
                    % pVals
                    toleg{end+1}=sprintf('%s vs %s',p_classes(iClass).trialName,p_classes(ic).trialName);
                    limval(1)=min([limval(1);comp(:)]);
                    limval(2)=max([limval(2);comp(:)]);
        
                end
            end
        end
    end
    title(totit{iChannel},'Interpreter','latex')
    % set(gca,'yscale','log')

    if ismember(newk(iChannel),selX)
        if length(time(inds))==1
            xlabel ('Model Comparison')
        else
            xlabel('time [ms]');
        end
    end
    if ismember(newk(iChannel),selY)
        % ylabel(par.ylabel,'interpreter','latex');
        ylabel('p (log scale)')
    end  
    
    yticks([0,0.01,0.05,0.1,0.5,1]);
    % legend(cc,toleg);
end
for iChannel=1:nChannels
    subplot(nRows_mod,nCols_mod,newk(iChannel))
    hold on; box on; grid on;
    set(gca,'yscale','log');
    ylim(limval);
    if length(time(inds))==1
        DX  = 0.5;
        xlim([1-DX,indcol+DX]);
    end
    % 
end

if legplot>0
  hLegend   = legend(cc,toleg);
  set(hLegend, 'position', get(hSub, 'position'));
  uistack(hLegend,'bottom');
end
p=[0,0,1500,600];
set(hfig,'Position',p);

if nargout>0
    h=hfig;
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end