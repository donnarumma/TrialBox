function distStruct = compute_pca_distances_clean(data_trials, par)

nTrials    = numel(data_trials);
InField    = par.InField;
dims       = par.wd;
centerIdx  = par.center;

distStruct = struct();

for iTr = 1:nTrials

    pcaMat = data_trials(iTr).(InField)(dims, :);
    time   = data_trials(iTr).(['time' InField]);

    refPoint = pcaMat(:, centerIdx);
    diffMat  = pcaMat - refPoint;
    dist     = sqrt(sum(diffMat.^2, 1));

    dt = diff(time);
    vel = [0, diff(dist) ./ dt];

    distStruct(iTr).distance  = dist;
    distStruct(iTr).velocity  = vel;
    distStruct(iTr).time      = time;
    distStruct(iTr).cond      = data_trials(iTr).Condition;
    distStruct(iTr).dir       = data_trials(iTr).Direction;
    distStruct(iTr).trialType = data_trials(iTr).trialType;
    distStruct(iTr).trialName = data_trials(iTr).trialName;
    distStruct(iTr).wd        = dims;

end
end