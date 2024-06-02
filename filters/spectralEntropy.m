function   [data_trials,out] = pEntropyCompute (data_trials,par)
% function [data_trials,out] = pEntropyCompute (data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
InField     = par.InField;
OutField    = par.OutField;
Inxfld      = 'time';
Outxfld     = 'time';
fact        = 8;

% data_pEntropy       = arrayfun(@(x)pentropy(x.(InField)',x.(['time' InField]))',data_trials,'UniformOutput',false);
for it=1:length(data_trials)    
    yData                               = data_trials(it).(InField);            % nVar x nTimes
    nVar                                = size(yData,1);
    xData                               = data_trials(it).([Inxfld InField]);   %    1 x nTimes
    nSamples                            = length(xData);               
    %% reverse engineering to get default tRes
    fs                                  = 1/mean(diff(xData));
    tRes                                = ceil(nSamples/fact)/fs;
    deltaT                              = tRes/4;
    nTimes                              = length((xData(1)+2*deltaT):(deltaT):(xData(end)-deltaT));
    %
    data_pEntropy                       = nan(nVar,nTimes);
    for iVar=1:nVar
        [data_pEntropy(iVar,:),tEst]    = pentropy(yData(iVar,:)',xData);
    end
    data_trials(it).(OutField)          = data_pEntropy;
    data_trials(it).([Outxfld OutField])= tEst';
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
