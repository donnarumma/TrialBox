function  [p_all,p_classes] = pSeparability(data_trials,par)
%function [p_all,p_classes] = pSeparability(data_trials,par)
execinfo        = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
trialTypes      = [data_trials.trialType];
trialNames      = {data_trials.trialName};
[classes,idx]   = unique(trialTypes);
names           = trialNames(idx);
xfld            = par.xfld;
InField         = par.InField;
OutField        = par.OutField;
time            = data_trials(1).([xfld InField]);
nChannels       = size(data_trials(1).(InField),1);
nTimes          = length(time);
nClasses        = length(classes);
nTrialsPerClass = length(data_trials(trialTypes==1));
X_data          = nan(nTrialsPerClass,nClasses);
p_data          = nan(nChannels,nTimes);
% p_classes       = cell(nClasses,1);
% comparisons-> nChannels x nTimes x nClasses -> stored
% pvalues    -> nChannels x nTimes            -> stored in 

for iTime=1:nTimes
    for iChannel=1:nChannels
        for ic=1:nClasses
            data_class=data_trials(trialTypes==classes(ic));
            for inc=1:nTrialsPerClass
                X_class         = data_class(inc).(InField);
                X_data(inc,ic)  = X_class(iChannel,iTime);
            end
        end
        [pval,~,st]             = anova1(X_data,[],'off');
        pval(pval==0)           = eps;
        comparison              = multcompare(st,'display','off');
        
        p_data(iChannel,iTime)  = pval;
        for iClass=1:nClasses
            p_classes(iClass).(OutField)(iChannel,iTime,:)=ones(1,1,nClasses);
            for iClassComp=1:nClasses
                if iClass == iClassComp
                    continue
                end
                ind = comparison(:,1)==iClass & comparison(:,2) == iClassComp;
                if sum(ind)==0
                    ind = comparison(:,2)==iClass & comparison(:,1) == iClassComp;
                end
                p_classes(iClass).(OutField)(iChannel,iTime,iClassComp)  = max(comparison(ind,end),eps);
                % COMPS(iClass,iChannel,iTime,iClassComp)               = comparison(ind,end);
            end
            p_classes(iClass).trialName           =names{iClass};
            p_classes(iClass).trialType           =classes(iClass);
            % Z_comps(iClass).(OutField2)         =Z_data;
            
            % Z_comps(iClass).([xfld OutField2])  =Z_comps(iClass).([xfld InField]);
        end
    end
end

p_all(1).(OutField)         = p_data;
vsTrialName=p_classes(1).trialName;
p_classes(1).([xfld OutField])=time;
for iClass=2:nClasses
    p_classes(iClass).([xfld OutField]) = time;
    vsTrialName = [vsTrialName ' vs ' p_classes(iClass).trialName];
end
p_all(1).([xfld OutField])   = data_trials.([xfld InField]);
p_all(1).trialName           = vsTrialName;
p_all(1).trialType           = 1;
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end

% %% Plot log(p)
% figure
% plot(time,log(p))

% [valmin,indmin]= min(log(pval));
% t_point = new_time(indmin);