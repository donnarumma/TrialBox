function   h=plot_EachDimVsTime(data_trials,par)
% function h=plot_EachDimVsTime(data_trials,par)
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
cmapslight      = par.cmapslight;

legplot         = par.legplot;
grayzone        = par.grayzone;
nCols           = par.nCols;
InField         = par.InField;
nComponents     = size(data_trials(1).(InField),1);
pos = get(hfig, 'position');
set(hfig, 'position', [pos(1) pos(2) 2*pos(3) pos(4)]);
xfld            = 'time';
timextkl        = data_trials.([xfld (InField)]);
nTimes          = size(timextkl,2);
nRows           = ceil(nComponents / nCols);
[~, idx]        = unique([data_trials.trialType]);
trialNames      = {data_trials.trialName};
names           = trialNames(idx);
if legplot==2
    hSub    = subplot(nRows+1, 1, 1);
    set(gca,'xticklabel','');set(gca,'yticklabel','');
    axis off;
    newk=indexesForSubplotsWithLegend2(nRows,nCols);
    nRows_mod=nRows+1;
    nCols_mod=nCols;
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
  
toleg=cell(0,0);
cc=[];
lw = 2.5;
for iTrial = 1:nTrials  %% NUMBER OF CLASSES
    col = cmaps(iTrial,:);
    col2= cmapslight(iTrial,:);
    for iComponent = 1:nComponents %% NUMBER OF COMPONENTS
    	if iTrial==1
        	ylcs{iComponent}=[+inf,-inf];
        end
        subplot(nRows_mod, nCols_mod, newk(iComponent));
        hold on;
        if ismember(iTrial,keep)    
            mtl = data_trials(iTrial).(InField)(iComponent,:);
            stl = data_trials(iTrial).([InField 'sigma'])(iComponent,:);
            if par.novariance
                stl=stl*0;
            end
            ylcs{iComponent}=[min([min(mtl-stl),ylcs{iComponent}(1)]),max([max(mtl+stl),ylcs{iComponent}(2)])];
            par.times       =timextkl;
            par.c1          =col;
            par.c2          =col2;
            par.c3          =col;
            par.fill        =1;
            par.hfig        =hfig;
            forP            =plotMeanAndVariancePlus(mtl,stl,par);
            plot(timextkl, mtl, 'linewidth', lw, 'color', col);      
            if iTrial==keep(end)
                try
                    ylmh=par.YLIM;                    
                catch
                    ylmh=ylcs{iComponent};      
                end
                h=diff(ylmh);
                dymh=h/100;
              
                if length(grayzone)==2
                    mg=grayzone(1);
                    sg=grayzone(2);
                else
                    mg=mean(grayzone);
                    sg=std(grayzone);
                end
                x=mg-sg;
                  
                y=ylmh(1);
              
                w=sg*2;
                try
                    mh1=rectangle('position',[x,y,w,h],'facecolor',0.9*ones(3,1),'linestyle','none');
                catch
                end
                mh2=plot([mg,mg],ylmh,'k--');
                uistack(mh2,'bottom');
                try
                    uistack(mh1,'bottom');
                catch
                end
                if ~par.novariance
                    ylim([ylmh(1)-dymh,ylmh(2)+dymh]);
                end
            end
            if iComponent==1 
                % fprintf('%s\n',names{iTrial});
                toleg{end+1}=names{iTrial};
                if nTimes>1
                    cc(end+1)=forP.hmean;
                else
                    cc(end+1)=bbb;
                end
            end  %% FIRST COMPONENTS LEGEND        
        end  %% KEEP VISIBLE
    end  %% END NUMBER OF COMPONENTS          
end  %% END NUMBER OF CLASSES

if legplot>0 
%   hLegend = legend(cc,toleg);
end    
selY=newk(1):nCols_mod:max(newk);

selX=newk(end-nCols+1):newk(end);

%% common labels
for iComponent = 1:nComponents %% NUMBER OF COMPONENTS
	h = subplot(nRows_mod, nCols_mod, newk(iComponent));

    box on; grid on;
    dx=(timextkl(2)-timextkl(1))/2;
    xl=[timextkl(1)-dx,timextkl(end)+dx];        
    xlim(xl);
    try
        strtl=par.titles{iComponent};
        title(strtl, 'interpreter', 'latex', 'fontsize', 16);
    catch
    end
    if ismember(newk(iComponent),selX) && nTimes>1
        xlabel('Time (ms)');
    end
    if ismember(newk(iComponent),selY)
        ylabel(par.ylabel,'interpreter','latex');
    end
end %% END NUMBER OF COMPONENTS
  
if legplot>0
  hLegend = legend(cc,toleg);
  set(hLegend, 'position', get(hSub, 'position'));
  uistack(hLegend,'bottom');
end
if nargout>0
    h=hfig;
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end