function [data_trials, pivotInfo] = splitLFPbyMonkey(data_trials, par)
% splitLFPbyMonkey
% Split LFP channels into monkey S and K and create final trial-wise signals.
%
% INPUT
%   data_trials : struct array containing trial data
%   par         : parameter struct with fields:
%                 - par.method
%                 - par.negThr
%
% OUTPUT
%   data_trials : updated struct array with LFP_S_final / LFP_K_final
%   pivotInfo   : diagnostic struct (mainly used for method = 'all')

    % -------------------------
    % Input validation
    % -------------------------
    if nargin < 2
        error('splitLFPbyMonkey requires data_trials and par as inputs.');
    end

    if ~isfield(par, 'method')
        error('Missing par.method.');
    end

    if ~isfield(par, 'negThr') || isempty(par.negThr)
        error('Missing par.negThr.');
    end

    numTrials = numel(data_trials);
    pivotInfo = struct();

    % -------------------------
    % Case 1: method = ''all''
    % -------------------------
    if strcmpi(par.method, 'all')

        % Split channels into S and K
        iLFP_S = data_trials(1).iLFP_S;
        iLFP_K = data_trials(1).iLFP_K;

        all_iLFP = [iLFP_S(:); iLFP_K(:)];
        idxS_mask = ismember(all_iLFP, iLFP_S);
        idxK_mask = ismember(all_iLFP, iLFP_K);

        nLFP_S = sum(idxS_mask);
        nLFP_K = sum(idxK_mask);

        for iTrial = 1:numTrials
            LFPmat = data_trials(iTrial).LFP;
            data_trials(iTrial).LFP_S = LFPmat(idxS_mask, :);
            data_trials(iTrial).LFP_K = LFPmat(idxK_mask, :);
            data_trials(iTrial).nLFP_S = nLFP_S;
            data_trials(iTrial).nLFP_K = nLFP_K;
        end

        % Global z-score + data-driven pivot + sign flip
        monkeys = {'S', 'K'};

        for m = 1:numel(monkeys)
            monkey = monkeys{m};
            fieldName = ['LFP_' monkey];
            nChannels = data_trials(1).(['nLFP_' monkey]);

            if nChannels == 0
                pivotInfo.(monkey).pivotIdx = [];
                pivotInfo.(monkey).flipSign = [];
                pivotInfo.(monkey).corrMat = [];
                pivotInfo.(monkey).centrality = [];
                continue;
            end

            % 1) Concatenate all session data for this monkey
            totalSamples = sum(arrayfun(@(x) size(x.(fieldName), 2), data_trials));
            all_data = zeros(nChannels, totalSamples);

            currentPtr = 1;
            for iT = 1:numTrials
                currentSize = size(data_trials(iT).(fieldName), 2);
                all_data(:, currentPtr:currentPtr + currentSize - 1) = data_trials(iT).(fieldName);
                currentPtr = currentPtr + currentSize;
            end

            % 2) Global z-score per channel
            mu_global = mean(all_data, 2);
            std_global = std(all_data, 0, 2);
            std_global(std_global == 0) = 1;

            all_data_z = (all_data - mu_global) ./ std_global;

            % 3) Correlation matrix, pivot selection, and sign flip
            if nChannels > 1
                R = corrcoef(all_data_z', 'Rows', 'pairwise');

                Rtmp = abs(R);
                Rtmp(1:nChannels+1:end) = NaN;
                centrality = median(Rtmp, 2, 'omitnan');

                [~, pivotIdx] = max(centrality);

                flipSign = ones(nChannels, 1);
                for iCh = 1:nChannels
                    rho = R(pivotIdx, iCh);
                    if ~isnan(rho) && rho < par.negThr
                        flipSign(iCh) = -1;
                    end
                end
                flipSign(pivotIdx) = 1;
            else
                R = 1;
                centrality = 1;
                pivotIdx = 1;
                flipSign = 1;
            end

            % Save diagnostic information
            pivotInfo.(monkey).pivotIdx = pivotIdx;
            pivotInfo.(monkey).flipSign = flipSign;
            pivotInfo.(monkey).corrMat = R;
            pivotInfo.(monkey).centrality = centrality;

            % 4) Apply z-score and sign flip trial by trial, then average channels
            for iT = 1:numTrials
                currentRaw = data_trials(iT).(fieldName);
                currentLFP_z = (currentRaw - mu_global) ./ std_global;

                if nChannels > 1
                    correctedLFP = currentLFP_z .* flipSign;
                    finalLFP = mean(correctedLFP, 1, 'omitnan');
                else
                    finalLFP = currentLFP_z;
                end

                data_trials(iT).(['LFP_' monkey '_final']) = finalLFP;
                data_trials(iT).(['pivotIdx_' monkey]) = pivotIdx;
            end
        end

    % -------------------------
    % Case 2: method ~= ''all''
    % -------------------------
    else
        monkeys = {'S', 'K'};

        for iTrial = 1:numTrials
            LFPmat = data_trials(iTrial).LFP;

            if size(LFPmat, 1) ~= numel(monkeys)
                error(['Expected LFP to be %dxn in non-all mode, found %dx%d in trial %d. ', ...
                       'Current method: %s'], ...
                       numel(monkeys), size(LFPmat,1), size(LFPmat,2), ...
                       iTrial, par.method);
            end

            for m = 1:numel(monkeys)
                monkey = monkeys{m};
                data_trials(iTrial).(['LFP_' monkey '_final']) = LFPmat(m, :);
            end
        end
    end

    % -------------------------
    % Add S/K direction labels
    % -------------------------
    data_trials = findSKdirection(data_trials);
end