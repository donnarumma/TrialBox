function   [EEG_Trials,out]=getFittsEEG(ind,par)
% function [EEG_Trials,out]=getFittsEEG(ind,par)
fprintf('Function: %s ',mfilename);
f_names = dir(par.data_dir);
row_del = find([f_names.isdir]==1); %Questo io ho dovuto inserirlo perché si prendeva i file nascosti
f_names(row_del)= [];

execinfo=par.exec;
if ~isempty(execinfo);t=tic; end

strinfo=sprintf('Loading Subject %g: %s%s%s',ind,f_names(ind).folder,filesep,f_names(ind).name);
fprintf('%s\n',strinfo)
subdata=load(strcat(f_names(ind).folder,filesep,f_names(ind).name));
EEG_Trials=subdata.trialdata_noartifacts;
EEG_Trials.filename=f_names(ind).name;if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end