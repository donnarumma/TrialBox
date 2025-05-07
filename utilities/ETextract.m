% function data = ETextract(data)

function data = ETextract(data)

for iTr = 1:length(data)
    try 
        indTonset = find([data(iTr).Event] == 2024);
    catch 
        indTonset = find([data(iTr).Event] == 2014);
    end
    try
        indMO_ET = find([data(iTr).Event] == 3000);
    catch
        indMO_ET = [];
    end
    
    data(iTr).ET = data(iTr).Event(indMO_ET,3)-data(iTr).Event(indTonset,3);
    data(iTr).timeET = data(iTr).timeJXYEXY - data(iTr).Event(indTonset,2);
    
end