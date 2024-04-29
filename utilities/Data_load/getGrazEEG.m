function   [EEG_Trials,out]=getGrazEEG(ind,par)
% function [EEG_Trials,out]=getGrazEEG(ind,par)
fprintf('Function: %s ',mfilename);
f_names = dir(par.data_dir);
f_names([f_names.isdir]==1)= [];

execinfo=par.exec;
if ~isempty(execinfo);t=tic; end

strinfo=sprintf('Loading Subject %g: %s%s%s',ind,f_names(ind).folder,filesep,f_names(ind).name);
fprintf('%s\n',strinfo)
subdata=load(strcat(f_names(ind).folder,filesep,f_names(ind).name));
EEG_Trials=subdata.data_row;
EEG_Trials.filename=f_names(ind).name;
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end