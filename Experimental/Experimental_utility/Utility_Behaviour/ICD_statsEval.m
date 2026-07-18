function ICD_out = ICD_statsEval(ICD_trials)

nTrials = length(ICD_trials);

ICD_out.trial = struct( ...
    'trialId', {}, ...
    'trialType', {}, ...
    'cond', {}, ...
    'dir', {}, ...
    'trialName', {}, ...
    'ICD', {}, ...
    'time_rel', {}, ...
    'time_win', {}, ...
    'XY1', {}, ...
    'XY2', {}, ...
    'JXY_window', {}, ...
    'ICD_mean', {}, ...
    'ICD_max', {}, ...
    'ICD_auc', {}, ...
    'ICD_timeMax', {}, ...
    'MT_fast_abs', {}, ...
    'MT_fast_rel', {}, ...
    'fastCode', {}, ...
    'alignCode', {}, ...
    'check_MT', {} );

k = 0;

for i = 1:nTrials
    
    if isempty(ICD_trials(i).ICD)
        continue;
    end
    
    ICD      = ICD_trials(i).ICD;
    time_rel = ICD_trials(i).time_rel;
    
    if isempty(ICD) || isempty(time_rel)
        continue;
    end
    
    k = k + 1;
    
    [ICDmax, idxMax] = max(ICD);
    
    ICD_out.trial(k).trialId      = ICD_trials(i).trialId;
    ICD_out.trial(k).trialType    = ICD_trials(i).trialType;
    ICD_out.trial(k).cond         = ICD_trials(i).cond;
    ICD_out.trial(k).dir          = ICD_trials(i).dir;
    ICD_out.trial(k).trialName    = ICD_trials(i).trialName;
    
    ICD_out.trial(k).ICD          = ICD_trials(i).ICD;
    ICD_out.trial(k).time_rel     = ICD_trials(i).time_rel;
    ICD_out.trial(k).time_win     = ICD_trials(i).time_win;
    ICD_out.trial(k).XY1          = ICD_trials(i).XY1;
    ICD_out.trial(k).XY2          = ICD_trials(i).XY2;
    ICD_out.trial(k).JXY_window   = ICD_trials(i).JXY_window;
    
    ICD_out.trial(k).ICD_mean     = mean(ICD, 'omitnan');
    ICD_out.trial(k).ICD_max      = ICDmax;
    ICD_out.trial(k).ICD_auc      = trapz(time_rel, ICD);
    ICD_out.trial(k).ICD_timeMax  = time_rel(idxMax);
    
    ICD_out.trial(k).MT_fast_abs  = ICD_trials(i).MT_fast_abs;
    ICD_out.trial(k).MT_fast_rel  = ICD_trials(i).MT_fast_rel;
    ICD_out.trial(k).fastCode     = ICD_trials(i).fastCode;
    ICD_out.trial(k).alignCode    = ICD_trials(i).alignCode;
    ICD_out.trial(k).check_MT     = ICD_trials(i).check_MT;
end

% Matrici cond x dir
ICD_out.meanMatrix    = meanTrialsToMatrix(ICD_out.trial, 'ICD_mean', 3, 8);
ICD_out.maxMatrix     = meanTrialsToMatrix(ICD_out.trial, 'ICD_max', 3, 8);
ICD_out.aucMatrix     = meanTrialsToMatrix(ICD_out.trial, 'ICD_auc', 3, 8);
ICD_out.timeMaxMatrix = meanTrialsToMatrix(ICD_out.trial, 'ICD_timeMax', 3, 8);

end
