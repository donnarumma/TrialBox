function [ICD_trials, MT_fast_trials] = ICD_trace(data_trials)

nTrials = length(data_trials);

ICD_trials(1,nTrials) = struct( ...
    'ICD', [], ...
    'time_rel', [], ...
    'time_win', [], ...
    'XY1', [], ...
    'XY2', [], ...
    'JXY_window', [], ...
    'MT_fast_abs', NaN, ...
    'MT_fast_rel', NaN, ...
    'fastCode', NaN, ...
    'alignCode', NaN, ...
    'check_MT', NaN, ...
    'trialId', NaN, ...
    'trialType', NaN, ...
    'cond', NaN, ...
    'dir', NaN, ...
    'trialName', '' );

MT_fast_trials = nan(nTrials,1);

tol = 1e-6;
minSamples = 20;

for i = 1:nTrials
    
    % Only Joint condition: DIR 1-8, ACT S, ACT K <=> trialType 17:24
    if data_trials(i).trialType < 17 || data_trials(i).trialType > 24
        continue;
    end
    
    Event = data_trials(i).Event;
    
    ICD_trials(i).trialType = data_trials(i).trialType;
    ICD_trials(i).trialId   = data_trials(i).trialId;
    ICD_trials(i).cond      = data_trials(i).Condition;
    ICD_trials(i).dir       = data_trials(i).Direction;
    
    if isfield(data_trials, 'trialName')
        ICD_trials(i).trialName = data_trials(i).trialName;
    end
    
    % Movement onset events
    idx2107 = find(Event(:,1) == 2107, 1, 'first');
    idx2207 = find(Event(:,1) == 2207, 1, 'first');
    
    if isempty(idx2107) || isempty(idx2207)
        warning('Trial %d: missing 2107 or 2207', i);
        continue;
    end
    
    time2107_abs = Event(idx2107, 2);
    time2207_abs = Event(idx2207, 2);
    
    if time2207_abs < time2107_abs
        fastIdx     = idx2207;
        fastCode    = 2207;
        fastTimeAbs = time2207_abs;
    else
        fastIdx     = idx2107;
        fastCode    = 2107;
        fastTimeAbs = time2107_abs;
    end
    
    % Alignment reference: prefer 2024, fallback 2014
    idx2024 = find(Event(:,1) == 2024, 1, 'first');
    idx2014 = find(Event(:,1) == 2014, 1, 'first');
    
    if ~isempty(idx2024)
        alignIdx     = idx2024;
        alignCode    = 2024;
        alignTimeAbs = Event(idx2024, 2);
    elseif ~isempty(idx2014)
        alignIdx     = idx2014;
        alignCode    = 2014;
        alignTimeAbs = Event(idx2014, 2);
    else
        warning('Trial %d: missing both 2024 and 2014', i);
        continue;
    end
    
    fastTimeRel = fastTimeAbs - alignTimeAbs;
    MT_fast_trials(i) = fastTimeAbs;
    
    % Consistency check using Event column 3
    fastTimeCol3     = Event(fastIdx, 3);
    alignTimeCol3    = Event(alignIdx, 3);
    fastTimeRelCheck = fastTimeCol3 - alignTimeCol3;
    
    isConsistent = abs(fastTimeRelCheck - fastTimeRel) < tol;
    
    if ~isConsistent
        warning('Trial %d: mismatch between Event col2 and col3 for fast MT', i);
    end
    
    time_data = data_trials(i).timeJXYEXY;
    JXYEXY    = data_trials(i).JXYEXY;
    
    % Window: 200 ms after faster movement onset
    idx_window = time_data >= fastTimeRel & time_data <= (fastTimeRel + 0.2);
    
    if sum(idx_window) < minSamples
        warning('Trial %d: ICD window too short', i);
        continue;
    end
    
    % Assuming rows = signals, columns = time
    XY1 = JXYEXY(1:2, idx_window);
    XY2 = JXYEXY(5:6, idx_window);
    
    ICD_values = sqrt(sum((XY1 - XY2).^2, 1));
    time_win   = time_data(idx_window);
    time_rel   = time_win - fastTimeRel;
    
    ICD_trials(i).ICD         = ICD_values;
    ICD_trials(i).time_rel    = time_rel;
    ICD_trials(i).time_win    = time_win;
    ICD_trials(i).XY1         = XY1;
    ICD_trials(i).XY2         = XY2;
    ICD_trials(i).JXY_window  = JXYEXY(:, idx_window);
    
    ICD_trials(i).MT_fast_abs = fastTimeAbs;
    ICD_trials(i).MT_fast_rel = fastTimeRel;
    ICD_trials(i).fastCode    = fastCode;
    ICD_trials(i).alignCode   = alignCode;
    ICD_trials(i).check_MT    = isConsistent;
end

end
