function dataOut = trimTrialsByTime(dataIn, params)
% TRIMTRIALSBYTIME Trim trial data and time vector to a specified time window.
%
%   dataOut = TRIMTRIALSBYTIME(dataIn, params) returns a copy of the input
%   struct array dataIn where both the data field specified by
%   params.dataField and its corresponding time vector ['time' dataField]
%   are cropped between params.tStart and params.tStop.
%
%   Inputs:
%     dataIn   - 1xN struct array, each element is a trial with fields:
%                  * timeX : time vector (e.g. 'timeSignal', 'timeLFP')
%                  * X      : data matrix [nChannels x nTime]
%                where X name is given by params.dataField.
%
%     params   - struct with fields:
%                  * dataField : char/string, name of the data field
%                                (e.g. 'Signal', 'LFP').
%                  * tStart    : scalar, start time of the window.
%                  * tStop     : scalar, stop time of the window.
%
%   Output:
%     dataOut  - struct array, same as dataIn but with:
%                  * dataField cropped in columns
%                  * corresponding timeX cropped consistently.

    % Extract parameters
    InField   = params.InField;
    tStart    = params.tStart;
    tStop     = params.tStop;

    nTrials = numel(dataIn);
    dataOut = dataIn;  

    % Build time field name and get time vector from first trial
    timeFieldName = ['time' InField];
    timeVector    = dataIn(1).(timeFieldName);

    % Find closest indices to tStart and tStop
    [~, idxStart] = min(abs(timeVector - tStart));
    [~, idxStop]  = min(abs(timeVector - tStop));

    % Ensure increasing indices
    if idxStart > idxStop
        tmp      = idxStart;
        idxStart = idxStop;
        idxStop  = tmp;
    end

    % Precompute cropped time vector
    croppedTime = timeVector(idxStart:idxStop);

    % Crop data and time for each trial
    for iTr = 1:nTrials
        trialData = dataIn(iTr).(InField);
        dataOut(iTr).(InField)   = trialData(:, idxStart:idxStop);
        dataOut(iTr).(timeFieldName) = croppedTime;
    end
end
