function data_trials = alignBattaglia(data_trials,par)

AlignEvent = par.AlignEvent;
InField = par.InField;
OutField = par.OutField;
xfld        = 'time';
xfld_old    = 'oldTime'; 

TimeLock = nan(1,length(data_trials));
for iTr = 1:length(data_trials)
    TimeLock(iTr) = data_trials(iTr).Event (data_trials(iTr).Event (:,1)==AlignEvent,2);
end

for iTr = 1:length(data_trials)
    data_trials(iTr).([xfld_old InField]) = data_trials(iTr).([xfld InField]);
    data_trials(iTr).([xfld OutField])  = data_trials(iTr).([xfld InField])  - TimeLock(iTr);
end
