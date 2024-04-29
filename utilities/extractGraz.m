% This is a function to extract Graz data for each subject
function [EEG_trial,fsample] = extractGraz(indsub,signal_name)

par.getGrazEEG.data_dir         = getDataGraz;
par.getGrazEEG.exec             = true;

EEG_raw                         = getGrazEEG(indsub,par.getGrazEEG);
fsample                         = EEG_raw.fsample; % Hertz

% FITTSarrangeTrials extract trials in array format in all time interval [-4,1]
par.GRAZarrangeTrials.InField  = signal_name; % name of the field in which to save data
par.GRAZarrangeTrials.it1      = -4;
par.GRAZarrangeTrials.it2      = 4;
par.GRAZarrangeTrials.fsample  = fsample;
par.GRAZarrangeTrials.exec     = true;
EEG_trial                      = GRAZarrangeTrials(EEG_raw.EEG_train,par.GRAZarrangeTrials);
