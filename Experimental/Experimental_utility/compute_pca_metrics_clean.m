function stats = compute_pca_metrics_clean(distStruct)

nTrials = numel(distStruct);

stats = struct();
stats.trial   = struct([]);
stats.condDir = struct([]);

allCond = nan(1, nTrials);
allDir  = nan(1, nTrials);

for i = 1:nTrials
    dist = distStruct(i).distance;
    vel  = distStruct(i).velocity;
    time = distStruct(i).time;

    postIdx = time >= 0;

    postTime = time(postIdx);
    postDist = dist(postIdx);
    postVel  = vel(postIdx);

    [MND, idxMND] = max(postDist);
    MNDT = postTime(idxMND);

    [dMND, idxdMND] = max(postVel);
    dMNDT = postTime(idxdMND);

    stats.trial(i).MND       = MND;
    stats.trial(i).MNDT      = MNDT;
    stats.trial(i).dMND      = dMND;
    stats.trial(i).dMNDT     = dMNDT;
    stats.trial(i).cond      = distStruct(i).cond;
    stats.trial(i).dir       = distStruct(i).dir;
    stats.trial(i).trialType = distStruct(i).trialType;
    stats.trial(i).trialName = distStruct(i).trialName;
    stats.trial(i).wd        = distStruct(i).wd;

    allCond(i) = distStruct(i).cond;
    allDir(i)  = distStruct(i).dir;
end

conds = unique(allCond);
dirs  = unique(allDir);

k = 1;
for c = conds
    for d = dirs
        idx = (allCond == c) & (allDir == d);

        if ~any(idx)
            continue
        end

        valsMND   = [stats.trial(idx).MND];
        valsMNDT  = [stats.trial(idx).MNDT];
        valsdMND  = [stats.trial(idx).dMND];
        valsdMNDT = [stats.trial(idx).dMNDT];

        stats.condDir(k).cond = c;
        stats.condDir(k).dir  = d;

        stats.condDir(k).MND   = valsMND;
        stats.condDir(k).MNDT  = valsMNDT;
        stats.condDir(k).dMND  = valsdMND;
        stats.condDir(k).dMNDT = valsdMNDT;

        k = k + 1;
    end
end

end