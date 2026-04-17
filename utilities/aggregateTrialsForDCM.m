function DCM_data = aggregateTrialsForDCM(data_trials, par)
% AGGREGATETRIALSFORDCM
% Aggregate single-trial data into Condition x Direction entries for DCM.
%
% INPUT
%   data_trials : struct array with one element per trial
%   par         : struct with at least:
%                 - par.session_name
%
% OUTPUT
%   DCM_data    : struct array with one element per Condition x Direction
%
% NOTES
%   - Old LFP is saved as startLFP.
%   - Old JXYEXY is saved as startJXYEXY.
%   - Old timeJXYEXY is saved as timestartJXYEXY.
%   - JXYEXY_LFP_fs is averaged across trials and saved as JXYEXY.
%   - Final tLFP and tJXYEXY are both equal to the final timeLFP.

    uniqueCond = unique([data_trials.Condition]);
    uniqueDir  = unique([data_trials.Direction]);

    DCM_data = struct();
    monkeys = {'S', 'K'};
    outIdx = 1;

    for iCond = 1:numel(uniqueCond)
        for iDir = 1:numel(uniqueDir)

            condValue = uniqueCond(iCond);
            dirValue  = uniqueDir(iDir);

            trialIdx = find([data_trials.Condition] == condValue & ...
                            [data_trials.Direction] == dirValue);

            if isempty(trialIdx)
                continue;
            end

            currentTrials = data_trials(trialIdx);
            nTrials = numel(currentTrials);

            % -------------------------------------------------------------
            % Metadata
            % -------------------------------------------------------------
            DCM_data(outIdx).Condition  = condValue;
            DCM_data(outIdx).Direction  = dirValue;
            DCM_data(outIdx).trialIds   = [currentTrials.trialId];
            DCM_data(outIdx).trialTypes = [currentTrials.trialType];
            DCM_data(outIdx).RT_S       = [currentTrials.RT_S];
            DCM_data(outIdx).RT_K       = [currentTrials.RT_K];
            DCM_data(outIdx).nTrials    = nTrials;
            DCM_data(outIdx).session    = par.session_name;

            if isfield(currentTrials, 'Chamber')
                DCM_data(outIdx).Chamber = currentTrials(1).Chamber;
            end
            if isfield(currentTrials, 'klabels')
                DCM_data(outIdx).klabels = currentTrials(1).klabels;
            end
            if isfield(currentTrials, 'nLFP_S')
                DCM_data(outIdx).nLFP_S = currentTrials(1).nLFP_S;
            end
            if isfield(currentTrials, 'nLFP_K')
                DCM_data(outIdx).nLFP_K = currentTrials(1).nLFP_K;
            end
            if isfield(currentTrials, 'selIndex')
                DCM_data(outIdx).selIndex = currentTrials(1).selIndex;
            end
            if isfield(currentTrials, 'iLFP_S')
                DCM_data(outIdx).iLFP_S = currentTrials(1).iLFP_S;
            end
            if isfield(currentTrials, 'iLFP_K')
                DCM_data(outIdx).iLFP_K = currentTrials(1).iLFP_K;
            end

            % -------------------------------------------------------------
            % Start / original signals
            % -------------------------------------------------------------
            if isfield(currentTrials, 'LFP')
                DCM_data(outIdx).startLFP = {currentTrials.LFP}.';
            else
                DCM_data(outIdx).startLFP = {};
            end

            if isfield(currentTrials, 'JXYEXY')
                DCM_data(outIdx).startJXYEXY = {currentTrials.JXYEXY}.';
            else
                DCM_data(outIdx).startJXYEXY = {};
            end

            if isfield(currentTrials, 'JXYEXY_LFP_fs')
                DCM_data(outIdx).startJXYEXY_LFP_fs = {currentTrials.JXYEXY_LFP_fs}.';
            else
                DCM_data(outIdx).startJXYEXY_LFP_fs = {};
            end
            
            if isfield(currentTrials, 'timeJXYEXY')
                DCM_data(outIdx).timestartJXYEXY = {currentTrials.timeJXYEXY}.';
            else
                DCM_data(outIdx).timestartJXYEXY = {};
            end

            % -------------------------------------------------------------
            % Keep existing fields unchanged
            % -------------------------------------------------------------
            if isfield(currentTrials, 'LFP_S_final')
                DCM_data(outIdx).OriginalLFP_S_final = {currentTrials.LFP_S_final}.';
            else
                DCM_data(outIdx).OriginalLFP_S_final = {};
            end

            if isfield(currentTrials, 'LFP_K_final')
                DCM_data(outIdx).OriginalLFP_K_final = {currentTrials.LFP_K_final}.';
            else
                DCM_data(outIdx).OriginalLFP_K_final = {};
            end

            if isfield(currentTrials, 'timeLFP')
                DCM_data(outIdx).OriginalTimeLFP_S_final = {currentTrials.timeLFP}.';
                DCM_data(outIdx).OriginalTimeLFP_K_final = {currentTrials.timeLFP}.';
            else
                DCM_data(outIdx).OriginalTimeLFP_S_final = {};
                DCM_data(outIdx).OriginalTimeLFP_K_final = {};
            end

            % -------------------------------------------------------------
            % Aggregate final LFP signals by monkey
            % -------------------------------------------------------------
            for iMonkey = 1:numel(monkeys)
                monkeyName = monkeys{iMonkey};
                signalField = ['LFP_' monkeyName '_final'];
                timeFieldOut = ['timeLFP_' monkeyName];
                signalFieldOut = ['LFP_' monkeyName];

                if ~isfield(currentTrials, signalField)
                    DCM_data(outIdx).(timeFieldOut) = [];
                    DCM_data(outIdx).(signalFieldOut) = [];
                    continue;
                end

                signalLengths = arrayfun(@(x) numel(x.(signalField)), currentTrials);
                [maxLen, longestTrialIdx] = max(signalLengths);

                if isfield(currentTrials, 'timeLFP')
                    DCM_data(outIdx).(timeFieldOut) = currentTrials(longestTrialIdx).timeLFP;
                else
                    DCM_data(outIdx).(timeFieldOut) = [];
                end

                paddedSignals = NaN(nTrials, maxLen);

                for iTrial = 1:nTrials
                    currentSignal = currentTrials(iTrial).(signalField);
                    currentLen = numel(currentSignal);
                    paddedSignals(iTrial, 1:currentLen) = currentSignal;
                end

                DCM_data(outIdx).(signalFieldOut) = mean(paddedSignals, 1, 'omitnan');
            end

            % -------------------------------------------------------------
            % Aggregate JXYEXY_LFP_fs across trials
            % Each row is averaged independently with NaN padding
            % -------------------------------------------------------------
            if isfield(currentTrials, 'JXYEXY_LFP_fs') && ~isempty(currentTrials(1).JXYEXY_LFP_fs)

                jxyLengths = arrayfun(@(x) size(x.JXYEXY_LFP_fs, 2), currentTrials);
                [maxLenJXY, ~] = max(jxyLengths);

                paddedRows = NaN(nTrials, 8, maxLenJXY);

                for iTrial = 1:nTrials
                    currentJXY = currentTrials(iTrial).JXYEXY_LFP_fs;

                    if isempty(currentJXY)
                        continue;
                    end

                    currentLen = size(currentJXY, 2);
                    nRowsCurrent = min(size(currentJXY, 1), 8);

                    paddedRows(iTrial, 1:nRowsCurrent, 1:currentLen) = ...
                        currentJXY(1:nRowsCurrent, 1:currentLen);
                end

                meanJXYEXY = squeeze(mean(paddedRows, 1, 'omitnan'));

                if isvector(meanJXYEXY)
                    meanJXYEXY = reshape(meanJXYEXY, 8, []);
                end

                DCM_data(outIdx).JXYEXY = meanJXYEXY;
            else
                DCM_data(outIdx).JXYEXY = [];
            end

            % -------------------------------------------------------------
            % Common time base
            % -------------------------------------------------------------
            if isfield(DCM_data(outIdx), 'timeLFP_S') && isfield(DCM_data(outIdx), 'timeLFP_K')
                if isequal(DCM_data(outIdx).timeLFP_S, DCM_data(outIdx).timeLFP_K)
                    DCM_data(outIdx).timeLFP = DCM_data(outIdx).timeLFP_S;
                else
                    if numel(DCM_data(outIdx).timeLFP_S) >= numel(DCM_data(outIdx).timeLFP_K)
                        DCM_data(outIdx).timeLFP = DCM_data(outIdx).timeLFP_S;
                    else
                        DCM_data(outIdx).timeLFP = DCM_data(outIdx).timeLFP_K;
                    end
                end
            else
                DCM_data(outIdx).timeLFP = [];
            end

            DCM_data(outIdx).tLFP = DCM_data(outIdx).timeLFP;
            DCM_data(outIdx).timeJXYEXY = DCM_data(outIdx).timeLFP;
            DCM_data(outIdx).tJXYEXY = DCM_data(outIdx).timeLFP;

            % -------------------------------------------------------------
            % Combined LFP matrix for DCM
            % -------------------------------------------------------------
            if isfield(DCM_data(outIdx), 'LFP_S') && isfield(DCM_data(outIdx), 'LFP_K') && ...
               ~isempty(DCM_data(outIdx).LFP_S) && ~isempty(DCM_data(outIdx).LFP_K)

                DCM_data(outIdx).LFP = [
                    DCM_data(outIdx).LFP_S;
                    DCM_data(outIdx).LFP_K
                ];
            else
                DCM_data(outIdx).LFP = [];
            end

            outIdx = outIdx + 1;
        end
    end
end