function   [data,out] = cspEncode(data,par)
% function [data,out] = cspEncode(data,par)
% csp executed per class
% inspired by csp algorithm https://ieeexplore.ieee.org/abstract/document/5332383

execinfo=par.exec;
if ~isempty(execinfo); t=tic; end

InField             = par.InField;
OutField            = par.OutField;
W                   = par.W;

data_project         = cat(4,data.(InField));

%% CSP Projection apply projection matrices on data
V_all                = CSPProjection(data_project,W);
data_CSP= cat(2,V_all{:});

%% organize data
xfld    = 'time';
for iTr=1:length(data)
    data(iTr).(OutField)            = data_CSP(iTr,:);
    data(iTr).([xfld OutField])=data(iTr).([xfld InField])(end);
end

%% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end
