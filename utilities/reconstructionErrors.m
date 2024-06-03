function   out = reconstructionErrors (data_trials,par)
% function out = reconstructionErrors (data_trials,par)
try par.exec; 
catch 
    par.exec=true; 
end
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

InFieldRec  = par.FieldReconstructed;
InFieldTar  = par.FieldTarget;

dimRec      = length(size(data_trials(1).(InFieldRec)));
dimTar      = length(size(data_trials(1).(InFieldTar)));

XdataRec    = cat(dimRec+1,data_trials.(InFieldRec));
XdataTar    = cat(dimTar+1,data_trials.(InFieldTar));
out         = errors(XdataRec,XdataTar); % output rec vs target rec
%% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end