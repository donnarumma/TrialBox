function data_trials = trimAlignedField(data_trials, par)
% trimAlignedField
% Trim an already aligned signal between tStart and min(tEnd, ET).
%
% INPUT
%   data_trials : struct array containing trial data
%   par         : struct with fields:
%                 - par.InField : signal field name, e.g. 'LFP', 'JXYEXY'
%                 - par.tStart  : start time for trimming
%                 - par.tEnd    : maximum end time for trimming
%
% OUTPUT
%   data_trials : updated struct array with trimmed signal and time fields

    inField   = par.InField;
    timeField = ['time' inField];
    tField    = ['t' inField];

    tStart = par.tStart;
    tEnd   = par.tEnd;

    for iTrial = 1:numel(data_trials)

        time_app = data_trials(iTrial).(timeField);
        data_app = data_trials(iTrial).(inField);

        idx_zero  = find(time_app >= 0, 1, 'first');
        idx_start = find(time_app >= tStart, 1, 'first');

        if isempty(idx_zero) || isempty(idx_start)
            continue
        end

        ET_absolute = data_trials(iTrial).ET;
        t_fin_effective = min(tEnd, ET_absolute);

        idx_end = find(time_app <= t_fin_effective, 1, 'last');

        if isempty(idx_end) || idx_start > idx_end
            continue
        end

        if isvector(data_app)
            data_trials(iTrial).(inField) = data_app(idx_start:idx_end);
        else
            data_trials(iTrial).(inField) = data_app(:, idx_start:idx_end);
        end

        data_trials(iTrial).(timeField) = time_app(idx_start:idx_end);
        % data_trials(iTrial).(tField)    = time_app(idx_start:idx_end);
    end
end