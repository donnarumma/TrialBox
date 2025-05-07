% function data = RTextract(data)

function data = RTextract(data)

for iTr = 1:length(data)
    try 
        indTonset = find([data(iTr).Event] == 2024);
    catch 
        indTonset = find([data(iTr).Event] == 2014);
    end
    try
        indMO_S = find([data(iTr).Event] == 2107);
    catch
        indMO_S = [];
    end
    try
        indMO_K = find([data(iTr).Event] == 2207);
    catch
        indMO_K = [];
    end

    if ~isempty(indMO_S) && isempty(indMO_K)
        data(iTr).RT_S = data(iTr).Event(indMO_S,3)-data(iTr).Event(indTonset,3);
        data(iTr).RT_K = NaN;
    elseif isempty(indMO_S) && ~isempty(indMO_K)
        data(iTr).RT_S = NaN;
        data(iTr).RT_K = data(iTr).Event(indMO_K,3)-data(iTr).Event(indTonset,3);
    elseif ~isempty(indMO_S) && ~isempty(indMO_K)
        data(iTr).RT_S = data(iTr).Event(indMO_S,3)-data(iTr).Event(indTonset,3);
        data(iTr).RT_K = data(iTr).Event(indMO_K,3)-data(iTr).Event(indTonset,3);
    end
end
