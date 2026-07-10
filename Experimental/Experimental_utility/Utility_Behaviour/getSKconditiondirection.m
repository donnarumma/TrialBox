function data = getSKconditiondirection(data)

for iTr=1:length(data)
    if contains(data(iTr).trialName,'D1') || ...
            contains(data(iTr).trialName,'DIR 1')
        data(iTr).Direction = 1;
    elseif contains(data(iTr).trialName,'D2')|| ...
            contains(data(iTr).trialName,'DIR 2')
        data(iTr).Direction = 2;
    elseif contains(data(iTr).trialName,'D3')|| ...
            contains(data(iTr).trialName,'DIR 3')
        data(iTr).Direction = 3;
    elseif contains(data(iTr).trialName,'D4')|| ...
            contains(data(iTr).trialName,'DIR 4')
        data(iTr).Direction = 4;
    elseif contains(data(iTr).trialName,'D5')|| ...
            contains(data(iTr).trialName,'DIR 5')
        data(iTr).Direction = 5;
    elseif contains(data(iTr).trialName,'D6')|| ...
            contains(data(iTr).trialName,'DIR 6')
        data(iTr).Direction = 6;
    elseif contains(data(iTr).trialName,'D7')|| ...
            contains(data(iTr).trialName,'DIR 7')
        data(iTr).Direction = 7;
    else
        data(iTr).Direction = 8;
    end

    if contains(data(iTr).trialName,'ACT S, OBS K')
       data(iTr).Condition = 1;
       data(iTr).ConditionName = 'Act S-Obs K';
    elseif contains(data(iTr).trialName,'OBS S, ACT K')
        data(iTr).Condition = 2;
        data(iTr).ConditionName = 'Obs S-Act K';
    else 
        data(iTr).Condition = 3;
        data(iTr).ConditionName = 'Act S-Act K';
    end       
end