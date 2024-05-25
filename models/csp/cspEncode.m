function   [data,out] = cspEncode(data,par)
% function [data,out] = cspEncode(data,par)
% csp executed per class
% inspired by csp algorithm https://ieeexplore.ieee.org/abstract/document/5332383

execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ', mfilename); end

InField             = par.InField;
OutField            = par.OutField;
xfld                = par.xfld;
W                   = par.W;

data_project         = cat(4,data.(InField));

% CSP Projection apply projection matrices on data
V_all                = cspProjection(data_project,W);
data_CSP= cat(2,V_all{:});

% organize data
for iTr=1:length(data)
    data(iTr).(OutField)            = data_CSP(iTr,:);
    data(iTr).([xfld OutField])=data(iTr).([xfld InField])(end);
end

% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
