function Group = averageBehaviorAcrossSessions(data_behaviour)
% AVERAGEBEHAVIORACROSSSESSIONS
%
% Compute population-level behavioural averages across SUA sessions.
%
% The function averages all scalar behavioural metrics stored in each
% session structure while ignoring NaN values.
%
% INPUT
% -------------------------------------------------------------------------
% data_behaviour : struct
%
%   Structure containing one field per session:
%
%       data_behaviour.SessionName
%
%   Each session must contain:
%
%       .S(1:24)
%           Behavioural metrics for monkey S
%
%       .K(1:24)
%           Behavioural metrics for monkey K
%
%   Optionally:
%
%       .ICD(1:24)
%           Inter-Cursor Distance metrics
%
%   Each element may contain scalar fields such as:
%
%       RT
%       PV
%       PVT
%       AE
%       AE_abs
%       ExT
%       EC
%       CMT
%       ICD_mean
%       ICD_max
%       ICD_auc
%
% OUTPUT
% -------------------------------------------------------------------------
% Group : struct
%
%   Group.S(1:24)
%       Mean behavioural metrics across sessions for monkey S.
%
%   Group.K(1:24)
%       Mean behavioural metrics across sessions for monkey K.
%
%   Group.ICD(1:24)
%       Mean ICD metrics across sessions.
%
% NOTES
% -------------------------------------------------------------------------
% - Only scalar numeric fields are averaged.
% - Missing values (NaN) are ignored.
% - Missing fields are skipped automatically.
%
% Author: Mirco Frosolone
% -------------------------------------------------------------------------

%% Validate input
if ~isstruct(data_behaviour)
    error('Input must be a struct containing session fields.');
end

sessNames = fieldnames(data_behaviour);

if isempty(sessNames)
    error('No sessions found in data_behaviour.');
end

nSess = numel(sessNames);

fprintf('\n');
fprintf('==============================================\n');
fprintf('Averaging behavioural metrics across sessions\n');
fprintf('Number of sessions: %d\n', nSess);
fprintf('==============================================\n\n');
%% Check session content
for s = 1:nSess

    sess = data_behaviour.(sessNames{s});

    if ~isstruct(sess)
        error('Session "%s" is not a valid structure.',sessNames{s});
    end

    if ~isfield(sess,'S') && ~isfield(sess,'K') && ~isfield(sess,'ICD')
        warning('Session "%s" contains no recognised behavioural fields.',...
            sessNames{s});
    end

end
%% Initialize output structure
Group = struct();
nEntries = 24; % 3 conditions x 8 directions
% Entry mapping:
%   1:8   -> Solo-S  (Act S / Obs K)
%   9:16  -> Solo-K  (Obs S / Act K)
%   17:24 -> Joint   (Act S / Act K)
Group.S = repmat(struct(),1,nEntries);
Group.K = repmat(struct(),1,nEntries);
hasICD = any(cellfun(@(x) ...
    isfield(data_behaviour.(x),'ICD'), sessNames));
if hasICD
    Group.ICD = repmat(struct(),1,nEntries);
end

%% Average metrics across condition-direction entries
for i = 1:nEntries
    S_list = cell(1,nSess);
    K_list = cell(1,nSess);

    nS = 0;
    nK = 0;
    if hasICD
        ICD_list = cell(1,nSess);
        nICD = 0;
    end
    for s = 1:nSess
        sess = data_behaviour.(sessNames{s});
        if isfield(sess,'S') && numel(sess.S) >= i
            nS = nS + 1;
            S_list{nS} = sess.S(i);
        end
        if isfield(sess,'K') && numel(sess.K) >= i
            nK = nK + 1;
            K_list{nK} = sess.K(i);
        end
        if hasICD && isfield(sess,'ICD') && numel(sess.ICD) >= i
            nICD = nICD + 1;
            ICD_list{nICD} = sess.ICD(i);
        end
    end

    S_list = S_list(1:nS);
    K_list = K_list(1:nK);
% ---- S ----
        if ~isempty(S_list)
            fnS = {};
            for s = 1:numel(S_list)
                fnS = [fnS; fieldnames(S_list{s})];
            end
            fnS = unique(fnS,'stable');
            for f = 1:numel(fnS)
                vals = nan(1,numel(S_list));
                nv = 0;
                for s = 1:numel(S_list)
                    if isfield(S_list{s}, fnS{f})
                        v = S_list{s}.(fnS{f});
                        if isnumeric(v) && isscalar(v)
                            nv = nv + 1;
                            vals(nv) = v;
                        end
                    end
                end
                vals = vals(1:nv);
                if ~isempty(vals)
                    Group.S(i).(fnS{f}) = mean(vals, 'omitnan');
                end
            end
        end
% ---- K ----
if ~isempty(K_list)
    fnK = {};
    for s = 1:numel(K_list)
        fnK = [fnK; fieldnames(K_list{s})];
    end
    fnK = unique(fnK,'stable');
            for f = 1:numel(fnK)
                vals = nan(1,numel(K_list));
                nv = 0;
                for s = 1:numel(K_list)
                    if isfield(K_list{s}, fnK{f})
                        v = K_list{s}.(fnK{f});
                        if isnumeric(v) && isscalar(v)
                            nv = nv + 1;
                            vals(nv) = v;
                        end
                    end
                end
                vals = vals(1:nv);
                if ~isempty(vals)
                    Group.K(i).(fnK{f}) = mean(vals, 'omitnan');
                end
            end
        end
% ---- ICD ----
    if hasICD
        ICD_list = ICD_list(1:nICD);
        if ~isempty(ICD_list)
            fnICD = {};
            for s = 1:numel(ICD_list)
                fnICD = [fnICD; fieldnames(ICD_list{s})];
            end
            fnICD = unique(fnICD,'stable');
            for f = 1:numel(fnICD)
                vals = nan(1,numel(ICD_list));
                nv = 0;
                for s = 1:numel(ICD_list)
                    if isfield(ICD_list{s},fnICD{f})
                        v = ICD_list{s}.(fnICD{f});
                        if isnumeric(v) && isscalar(v)
                            nv = nv + 1;
                            vals(nv) = v;
                        end
                    end
                end
                vals = vals(1:nv);
                if ~isempty(vals)
                    Group.ICD(i).(fnICD{f}) = mean(vals,'omitnan');
                end
            end
        end
    end
end
fprintf('\n');
fprintf('==============================================\n');
fprintf('Behavioural averaging completed\n');
fprintf('Sessions processed : %d\n', nSess);
fprintf('Entries averaged   : %d\n', nEntries);
fprintf('Generated fields:\n');

if isfield(Group,'S')
    fprintf('   - Group.S\n');
end

if isfield(Group,'K')
    fprintf('   - Group.K\n');
end

if isfield(Group,'ICD')
    fprintf('   - Group.ICD\n');
end

fprintf('==============================================\n\n');
end