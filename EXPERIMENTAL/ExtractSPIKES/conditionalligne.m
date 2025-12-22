function data_trials = conditionalligne(data_trials)

for iTr=1:length(data_trials)
    if data_trials(iTr).trialTypeCond == 1
        data_trials(iTr).trialTypeCond = 3;
    elseif data_trials(iTr).trialTypeCond == 2
        data_trials(iTr).trialTypeCond = 1;
    elseif data_trials(iTr).trialTypeCond == 3
        data_trials(iTr).trialTypeCond = 2;
    end
end
