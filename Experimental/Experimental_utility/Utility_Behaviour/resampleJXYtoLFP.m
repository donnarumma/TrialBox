function data_trials = resampleJXYtoLFP(data_trials)
% resampleJXYtoLFP
% Resample JXYEXY (2 channels) to LFP time axis trial by trial.

    for iTrial = 1:numel(data_trials)
        
        % Skip if no data
        if isempty(data_trials(iTrial).LFP) || isempty(data_trials(iTrial).JXYEXY)
            continue;
        end
        
        time_jxy = data_trials(iTrial).timeJXYEXY;
        jxy_data = data_trials(iTrial).JXYEXY;
        time_lfp = data_trials(iTrial).timeLFP;
        
        % Interpolate each JXYEXY channel separately
        jxy_resampled = zeros(2, length(time_lfp));
        for ch = 1:size(jxy_data,1)
            jxy_resampled(ch, :) = interp1(time_jxy, jxy_data(ch, :), ...
                                         time_lfp, 'linear', 'extrap');
        end
        
        data_trials(iTrial).JXYEXY_LFP_fs = jxy_resampled;
    end
end