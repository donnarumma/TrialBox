
function data = findSKdirection(data)

for iTr=1:length(data)
    if ismember(data(iTr).trialType,1:8)
       data(iTr).ConditionName = 'Act S-Obs K';
       data(iTr).Condition = 1;
    elseif ismember(data(iTr).trialType,9:16)
        data(iTr).ConditionName = 'Obs S-Act K';
        data(iTr).Condition = 2;
    else 
        data(iTr).ConditionName = 'Act S-Act K';
        data(iTr).Condition = 3;
    end
    position = find(data(iTr).klabels ==1);
    if ismember(position,[1,9,17])
        data(iTr).Direction = 1;
    elseif ismember(position,[2,10,18])
        data(iTr).Direction = 2;
    elseif ismember(position,[3,11,19])
        data(iTr).Direction = 3;
    elseif ismember(position,[4,12,20])
        data(iTr).Direction = 4;
    elseif ismember(position,[5,13,21])
        data(iTr).Direction = 5;
    elseif ismember(position,[6,14,22])
        data(iTr).Direction = 6;
    elseif ismember(position,[7,15,23])
        data(iTr).Direction = 7;
    elseif ismember(position,[8,16,24])
        data(iTr).Direction = 8;
    end
end