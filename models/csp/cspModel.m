function   [data,out] = cspModel(data,par)
% function Out = cspModel(data,par)
% csp executed per class
% inspired by csp algorithm https://ieeexplore.ieee.org/abstract/document/5332383
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s', mfilename); end

numcomponents   = par.m;
InField         = par.InField;

data_4d         = cat(4,data.(InField));
labs            = [data.trialType]';

% CSP all data
out.W           = cspDictionary(data_4d,labs, numcomponents);

% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf(' | Time Elapsed: %.2f s\n',out.exectime); end
