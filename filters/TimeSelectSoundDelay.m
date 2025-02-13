function   [data_trials, out] = TimeSelectSoundDelay(data_trials,par)
% function [data_trials, out] = TimeSelectSoundDelay(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ', mfilename); end

dt          = par.dt;
InField     = par.InField;
OutField    = par.OutField;
xfld        = 'time';

%% time selection and interpolation: to put in a separate function
t1_delay = par.t1;
t2_delay = par.t2;
nTrials             = length(data_trials);
for iTrial = 1:nTrials
    t1 = t1_delay(iTrial);

    if par.deltaT
        t2 = t2_delay(iTrial) + t1 + par.deltaT;
    else
        t2 = t2_delay(iTrial);
    end

    [~,itstart]  = min(abs(data_trials(iTrial).([xfld InField]) - t1));
    [~,itstop]   = min(abs(data_trials(iTrial).([xfld InField]) - t2));
    data_trials(iTrial).(OutField)         = data_trials(iTrial).(InField)(:,itstart:itstop-1);
    data_trials(iTrial).([xfld OutField])  = data_trials(iTrial).([xfld InField])(itstart:itstop-1);
    
    if dt>0
        n_var                                   = size(data_trials(iTrial).(InField),1);
        X_inter=nan(n_var,length(itstart:dt:itstop));
        for i_var=1:n_var
            X_inter(i_var,:)=interp1(itstart:itstop,data_trials(iTrial).(OutField)(i_var,:),itstart:dt:itstop,'spline','extrap');
        end
        data_trials(iTrial).(OutField)         = X_inter;
        data_trials(iTrial).([xfld OutField])  = interp1(itstart:itstop,data_trials(iTrial).([xfld OutField]),itstart:dt:itstop,'spline','extrap');
    end
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Time Elapsed: %.2f s\n',out.exectime); end