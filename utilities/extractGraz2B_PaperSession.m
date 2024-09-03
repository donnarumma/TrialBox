% This is a function to extract Graz data for each subject
function [EEG_trial,fsample] = extractGraz2B_PaperSession(indsub,par)


InField                         = par.InField ;
signal_name                     = par.signal_name;

par.getGrazEEG.data_dir         = getDataGraz2B_PaperSession;
par.getGrazEEG.exec             = true;

EEG_raw                         = getGrazEEG(indsub,par.getGrazEEG);
fsample                         = EEG_raw.fsample; % Hertz

% FITTSarrangeTrials extract trials in array format in all time interval [-4,1]
par.GRAZarrangeTrials2B.InField  = signal_name; % name of the field in which to save data
par.GRAZarrangeTrials2B.it1      = -5;
par.GRAZarrangeTrials2B.it2      = 4.5;
par.GRAZarrangeTrials2B.fsample  = fsample;
par.GRAZarrangeTrials2B.exec     = true;
EEG_trial                      = GRAZarrangeTrials2B(EEG_raw.(strcat('EEG_',InField)),par.GRAZarrangeTrials2B);
