function data_trials = alignBattaglia(data_trials,par)

AlignEvent = par.AlignEvent;
LockEvent  = par.LockEvent;
InField    = par.InField;
OutField   = par.OutField;

xfld     = 'time';
xfld_old = 'originalTime';

for iTr = 1:length(data_trials)

    events = data_trials(iTr).Event;

    idxAlign = find(events(:,1) == AlignEvent, 1);
    idxLock  = find(events(:,1) == LockEvent, 1);

    if isempty(idxAlign)
        warning('AlignEvent non trovato nel trial %d', iTr);
        continue
    end

    tAlign = events(idxAlign,2);

    data_trials(iTr).([xfld_old InField]) = data_trials(iTr).([xfld InField]);
    data_trials(iTr).([xfld OutField])    = data_trials(iTr).([xfld InField]) - tAlign;

    if ~isempty(idxLock) && isfield(data_trials(iTr),'ET') && ~isempty(data_trials(iTr).ET)
        tLock = events(idxLock,2);
        delta = tLock - tAlign;
        data_trials(iTr).ET = data_trials(iTr).ET + delta;
    end
end