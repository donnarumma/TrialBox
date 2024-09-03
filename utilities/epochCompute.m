function   [data_out, out] = epochCompute(data_trials,params)
% function [data, out] = epochCompute(data,params)

execinfo=params.exec;
if ~isempty(execinfo); t=tic; end

InField           = params.InField;
fsample           = params.fample;
t_epoch           = params.t_epoch;
lengEpoch         = floor(t_epoch*fsample); 
overlap_percent   = params.overlap_percent;

data_out = data_trials;
num_chan = size(data_trials(1).(InField),1);
num_filt = size(data_trials(1).(InField),3);
time_intervals = struct();
sizetime = NaN(size(data_trials,1),1);
for nTrials = 1:size(data_trials,1)
    [n_epoch,time_intervals(nTrials).T,bin_intervals] = computeEpochIntervals(size(data_trials(nTrials).(InField),2)/fsample,t_epoch,overlap_percent,fsample);
    sizetime(nTrials,1) = size(time_intervals(nTrials).T,1);
    eeg_trial = data_trials(nTrials).(InField);
    EEG_final = NaN(n_epoch,num_chan,lengEpoch,num_filt);
    for nEp = 1:n_epoch
        id_start = bin_intervals(nEp,1);
        id_stop = bin_intervals(nEp,2);
        EEG_final(nEp,:,:,:)= eeg_trial(:,id_start:id_stop,:);
    end
    data_out(nTrials).(InField) = EEG_final;
end
out.time_intervals = time_intervals(max(sizetime)).T;
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Time Elapsed: %.2f s\n',out.exectime); end