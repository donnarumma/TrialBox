function Group = averageBehaviorAcrossSessions(data_behaviour)
% averageBehaviorAcrossSessions
% Input:
%   data_behaviour.session1, data_behaviour.session2, ...
%   each session is a struct with fields:
%     .S(1:24).<fields>
%     .K(1:24).<fields>
%
% Output:
%   Group.S(1:24).<fields> and Group.K(1:24).<fields>
%   as means across sessions (omit NaNs)

    if ~isstruct(data_behaviour)
        error('Input must be a struct containing session fields.');
    end

    sessNames = fieldnames(data_behaviour);
    nSess = numel(sessNames);
    if nSess == 0
        error('No sessions found in data_behaviour.');
    end

    Group = struct();
    Group.S = repmat(struct(),1,24);
    Group.K = repmat(struct(),1,24);

    for i = 1:24
        S_list = cell(1,nSess);
        K_list = cell(1,nSess);
        nS = 0;
        nK = 0;

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
        end

        S_list = S_list(1:nS);
        K_list = K_list(1:nK);

        if ~isempty(S_list)
            fnS = fieldnames(S_list{1});
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

        if ~isempty(K_list)
            fnK = fieldnames(K_list{1});
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
    end
end