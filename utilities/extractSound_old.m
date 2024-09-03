% This is a function to extract Sound data for each subject
function [EEG_trial,fsample] = extractSound_old(indsub,par)

signal_name                      = par.signal_name;
it_end                           = par.it_end; 

par.getSoundEEG.data_dir         = getDataSound_old;
par.getSoundEEG.exec             = true;

EEG_raw                          = getSoundEEG_old(indsub,par.getSoundEEG);
fsample                          = EEG_raw.fsample;

% FITTSarrangeTrials extract trials in array format in all time interval [-4,1]
par.SOUNDarrangeTrials.it_start  = 1/fsample;
par.SOUNDarrangeTrials.it_end    = it_end;
par.SOUNDarrangeTrials.InField   = signal_name; % name of the field in which to save data
par.SOUNDarrangeTrials.OutField  = signal_name; % name of the field in which to save data
par.SOUNDarrangeTrials.ch_number = EEG_raw.ch_number;
par.SOUNDarrangeTrials.exec      = true;

if par.multiEpoch == true
    EEG_trial                        = SOUNDarrangeTrialsMultiEpoch_OLD(EEG_raw,par.SOUNDarrangeTrials);
else 
    EEG_trial                        = SOUNDarrangeTrials(EEG_raw,par.SOUNDarrangeTrials);
end