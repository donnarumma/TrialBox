% This is a function to extract Sound data for each subject
function [EEG_trial,fsample] = extractSound(indsub,par)

signal_name                      = par.signal_name;

par.getSoundEEG.data_dir         = getDataSound;

par.getSoundEEG.exec             = true;

EEG_raw                          = getSoundEEG(indsub,par.getSoundEEG);
fsample                          = EEG_raw.fs;

% FITTSarrangeTrials extract trials in array format in all time interval [-6,9]
par.SOUNDarrangeTrials.InField   = signal_name; % name of the field in which to save data
par.SOUNDarrangeTrials.OutField  = signal_name; % name of the field in which to save data
par.SOUNDarrangeTrials.exec      = true;

if par.multiEpoch == true
    EEG_trial                        = SOUNDarrangeTrialsMultiEpoch(EEG_raw,par.SOUNDarrangeTrials);
else 
    EEG_trial                        = SOUNDarrangeTrials(EEG_raw,par.SOUNDarrangeTrials);
end