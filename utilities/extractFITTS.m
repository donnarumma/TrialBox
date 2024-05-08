% This is a function to extract FITTS data for each subject
function [EEG_trials,fsample] = extractFITTS(indsub,par)

signal_name                     = par.signal_name;
directoryInfield                = par.dir;

if strcmp(directoryInfield,'epogo')
    par.getFittsEEG.data_dir         = getDataFITTSepogo;
else
    par.getFittsEEG.data_dir         = getDataFITTSepotarget;
end
par.getFittsEEG.exec                 = true;

EEG_raw                              = getFittsEEG(indsub,par.getFittsEEG);
fsample                              = 1000; % Hertz

% FITTSarrangeTrials extract trials in array format in all time interval [-4,1]
t_TargON                             = find(round(EEG_raw.time{1},3)==0);
it1                                  = 1;              % first bin
it2                                  = t_TargON+1000-1;    % 1000 time bins after go

par.FITTSarrangeTrials.it_start      = it1;
par.FITTSarrangeTrials.it_end        = it2;
if strcmp(directoryInfield,'epogo')
    par.FITTSarrangeTrials.event         = 'GO';
else
    par.FITTSarrangeTrials.event         = 'TargetOFF';
end
par.FITTSarrangeTrials.InField       = signal_name; % name of the field in which to save data
par.FITTSarrangeTrials.ch_number     = 59;
par.FITTSarrangeTrials.exec          = true;
EEG_trials                           = FITTSarrangeTrials(EEG_raw,par.FITTSarrangeTrials);
