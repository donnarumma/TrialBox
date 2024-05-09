function   [EEG_Trials,out]=getSoundEEG(ind,par)
% function [EEG_Trials,out]=getSoundEEG(ind,par)
fprintf('Function: %s ',mfilename);
f_names = dir(par.data_dir);
f_names([f_names.isdir]==1)= [];

execinfo=par.exec;
if ~isempty(execinfo);t=tic; end

strinfo=sprintf('Loading Subject %g: %s%s%s',ind,f_names(ind).folder,filesep,f_names(ind).name);
fprintf('%s\n',strinfo)
subdata=load(strcat(f_names(ind).folder,filesep,f_names(ind).name));
EEG_Trials.filename         = f_names(ind).name;
EEG_Trials.eeg_sound        = subdata.EEG.X_sound;
EEG_Trials.eeg_text         = subdata.EEG.X_text;
EEG_Trials.eeg_baseline     = subdata.EEG.X_baseline;
EEG_Trials.ch_number        = subdata.EEG.nbchan;
EEG_Trials.fsample          = subdata.EEG.fs;
EEG_Trials.class_indices    = subdata.EEG.class_indices;
EEG_Trials.y                = subdata.EEG.y;
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end