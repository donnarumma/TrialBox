function   h=plot_EachDimBar(data_trials,par)
% function h=plot_EachDimBar(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig    =figure;
    par.hfig=hfig;
else
    hfig    =par.hfig;
end  
nTrials         = length(data_trials);
keep            = par.keep;
cmaps           = par.cmaps;
legplot         = par.legplot;
grayzone        = par.grayzone;
nCols           = par.nCols;
InField         = par.InField;

dp              = 1;
if isempty(par.addbar)
    cumulative  = false;
else
    cumulative  = true;
    if nCols>1
        dp          = 1.5;
    end
end
evaltime        = par.evaltime;

[~, idx]        = unique([data_trials.trialType]);
trialNames      = {data_trials.trialName};
names           = trialNames(idx);
nComponents     = size(data_trials(1).(InField),1);
pos             = get(hfig, 'position');
set(hfig, 'position', [pos(1) pos(2) 2*pos(3)*dp pos(4)]);
xfld            = 'time';
% 
nRows           = ceil(nComponents / nCols);
if nRows==1 && nCols==1
    nRows_mod   = nRows;
    nCols_mod   = nCols;
    newk        = 1;
elseif legplot==2
    hSub        = subplot(nRows+1, 1, 1);
    set(gca,'xticklabel','');set(gca,'yticklabel','');
    axis off;
    newk        = indexesForSubplotsWithLegend2(nRows,nCols);
    nRows_mod   = nRows+1;
    nCols_mod   = nCols;
elseif legplot==1
    hSub        = subplot(1, nCols+1, 1);
    set(gca,'xticklabel','');set(gca,'yticklabel','');
    axis off;
    newk        = indexesForSubplotsWithLegend(nRows,nCols);
    nRows_mod   = nRows;
    nCols_mod   = nCols+1;
end
  
toleg=cell(0,0);
cc=[];

% XXX=cell(nComponents,1);  
for iTrial = 1:nTrials  %% NUMBER OF CLASSES
    col = cmaps(iTrial,:);
    for iComponent = 1:nComponents %% NUMBER OF COMPONENTS
    	if iTrial==1
        	ylcs{iComponent}=[+inf,-inf];
        end
        subplot(nRows_mod, nCols_mod, newk(iComponent));
        hold on;
     
        if ismember(iTrial,keep) 
            timextkl        = data_trials.([xfld (InField)]);
            [~,indextime]   = min(abs(timextkl-evaltime));
            mtl = data_trials(iTrial).(InField)(iComponent,indextime);
            stl = data_trials(iTrial).([InField 'sigma'])(iComponent,indextime);
            if par.novariance
                stl=stl*0;
            end
            ylcs{iComponent}=[min([min(mtl-stl),ylcs{iComponent}(1)]),max([max(mtl+stl),ylcs{iComponent}(2)])];

            %% BAR CASE 
            bbb=bar(iTrial,mtl,'facecolor',col);
            errorbar(iTrial,mtl,stl,'ko','linewidth',2)
        
            if iTrial==keep(end)
                try
                    ylmh=par.YLIM;                    
                catch
                    ylmh=ylcs{iComponent};      
                end
                h=diff(ylmh);
                dymh=h/100;
              
                if length(grayzone)==2
                    mg=grayzone(1);%
                    sg=grayzone(2);%
                else
                    mg=mean(grayzone);
                    sg=std(grayzone);
                end
                x=mg-sg;
                  
                y=ylmh(1);
                if ~par.novariance
                    ylim([ylmh(1)-dymh,ylmh(2)+dymh]);
                end
            end
            if iComponent==1 
                % fprintf('%s\n',names{iTrial});
                toleg{end+1}=names{iTrial};
                cc(end+1)=bbb;
            end  %% FIRST COMPONENTS LEGEND        
        end  %% KEEP VISIBLE
    end  %% END NUMBER OF COMPONENTS          
end  %% END NUMBER OF CLASSES

selY=newk(1):nCols_mod:max(newk);

selX=newk(end-nCols+1):newk(end);

%% common labels
for iComponent = 1:nComponents %% NUMBER OF COMPONENTS
	h = subplot(nRows_mod, nCols_mod, newk(iComponent));

    box on; grid on;
    %% BAR CASE
    xl=[1-0.5,nTrials+0.5];
    if cumulative
        timextkl        = par.addbar(1).([xfld (InField)]);
        [~,indextime]   = min(abs(timextkl-evaltime));
        mtl             = par.addbar(1).(InField)(iComponent,indextime);
        stl             = par.addbar(1).([InField, 'sigma'])(iComponent,indextime);
        if par.novariance
            stl=stl*0;
        end
        errorbar(nTrials+2,mtl,stl,'ko','linewidth',2);
        bbb     =bar(nTrials+2,mtl,'facecolor',0.4*[1,1,1]);
        errorbar(nTrials+2,mtl,stl,'ko','linewidth',2);
        xl      =[1-0.5,nTrials+2+0.5];
        
        plot([nTrials+1,nTrials+1],[0,1],'k');
        if iComponent==1
            toleg{end+1}='CUMULATIVE';
            cc(end+1)=bbb;
        end
    end
    xlim(xl);
    try
        strtl=par.titles{iComponent};
        if cumulative
            % strtl=sprintf('%s - Tot\\%% $$%.2g',strtl,mtl);
            strtl=sprintf('%s - tot $$%.2g',strtl,mtl);
            if ~par.novariance
                % strtl=sprintf('%s\\pm%g$$',strtl,round(stl*100));
                strtl=sprintf('%s\\pm%.2g$$',strtl,stl);
            else
                strtl=sprintf('%s$$',strtl);
            end
        end
        title(strtl, 'interpreter', 'latex', 'fontsize', 16);
    catch
    end
    if ismember(newk(iComponent),selY)
        ylabel(par.ylabel,'interpreter','latex');
    end
    if par.chanceline 
        set(gca,'Xtick',[]);
        plot(xl,[1/nTrials,1/nTrials],'k--','LineWidth',1)
    end
end %% END NUMBER OF COMPONENTS
  
if legplot>0
    hLegend = legend(cc,toleg);
    if nRows_mod>1 || nCols_mod>1
        set(hLegend, 'position', get(hSub, 'position'));
        uistack(hLegend,'bottom');
    else
        if     legplot==1
            set(hLegend,'Location','WestOutside');
        elseif legplot==2
            set(hLegend,'Location','NorthOutside');
        end
    end
end
if nargout>0
    h=hfig;
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end