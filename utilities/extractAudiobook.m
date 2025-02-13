% This is a function to extract Audiobook data for each subject
function [EEG_trial,fsample] = extractAudiobook(indsub,par)

signal_name                      = par.signal_name;

par.getAudiobookEEG.data_dir         = getDataAudiobook;

par.getAudiobookEEG.exec             = true;

EEG_raw                          = getAudiobookEEG(indsub,par.getAudiobookEEG);
fsample                          = EEG_raw.fs;

% FITTSarrangeTrials extract trials in array format in all time interval [-6,9]
par.AUDIOBOOKarrangeTrials.InField   = signal_name; % name of the field in which to save data
par.AUDIOBOOKarrangeTrials.OutField  = signal_name; % name of the field in which to save data
par.AUDIOBOOKarrangeTrials.exec      = true;


if par.multiEpoch == true
    EEG_trial                        = AUDIOBOOKarrangeTrialsMultiEpoch(EEG_raw,par.AUDIOBOOKarrangeTrials);
else
    par.AUDIOBOOKarrangeTrials.it_start  = par.it_start;
    par.AUDIOBOOKarrangeTrials.it_end    = par.it_end;
    par.AUDIOBOOKarrangeTrials.fs        = fsample;

    EEG_trial                        = AUDIOBOOKarrangeTrials(EEG_raw,par.AUDIOBOOKarrangeTrials);
end