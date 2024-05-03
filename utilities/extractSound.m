% This is a function to extract Sound data for each subject
function [EEG_trial,fsample] = extractSound(indsub,par)

signal_name                      = par.signal_name;

par.getSoundEEG.data_dir         = getDataSound;
par.getSoundEEG.exec             = true;

EEG_raw                          = getSoundEEG(indsub,par.getSoundEEG);
fsample                          = EEG_raw.fsample;

% FITTSarrangeTrials extract trials in array format in all time interval [-4,1]
par.SOUNDarrangeTrials.it_in     = 1/fsample;
par.SOUNDarrangeTrials.it_sum    = 2.5;
par.SOUNDarrangeTrials.InField   = signal_name; % name of the field in which to save data
par.SOUNDarrangeTrials.OutField  = signal_name; % name of the field in which to save data
par.SOUNDarrangeTrials.ch_number = EEG_raw.ch_number;
par.SOUNDarrangeTrials.exec      = true;

EEG_trial                        = SOUNDarrangeTrials(EEG_raw,par.SOUNDarrangeTrials);