function   [data_trials,out] = instantaneousFrequency (data_trials,par)
% function [data_trials,out] = instantaneousFrequency (data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
InField     = par.InField;
OutField    = par.OutField;
Inxfld      = 'time';
Outxfld     = 'time';                   
fact        = 8;
% data_instfreq       = cellfun(@(x,t)instfreq(x,t)',data,time,'UniformOutput',false); 
% data_instfreq       = arrayfun(@(x)instfreq(x.(InField)',x.([Inxfld InField]))',data_trials,'UniformOutput',false);
for it=1:length(data_trials)    
    yData                               = data_trials(it).(InField);            % nVar x nTimes
    nVar                                = size(yData,1);
    xData                               = data_trials(it).([Inxfld InField]);   %    1 x nTimes
    nSamples                            = length(xData);               
    %%%%
    % data_instfreq                       = instfreq(yData',xData)';
    %% reverse engineering to get default tRes
    fs                                  = 1/mean(diff(xData));
    tRes                                = ceil(nSamples/fact)/fs;
    deltaT                              = tRes/4;
    nTimes                              = length((xData(1)+2*deltaT):(deltaT):(xData(end)-deltaT));
    %
    data_instfreq                       = nan(nVar,nTimes);
    for iVar=1:nVar
        [p,f,tEst]                      = pspectrum(yData(iVar,:),xData,'spectrogram','TimeResolution',tRes);
        data_instfreq(iVar,:)           = instfreq(p,f,tEst-tEst(1));
    end
    data_trials(it).(OutField)           = data_instfreq;
   
    % timein                               = data_trials(it).([Inxfld InField]);
    % T                                    = size(data_trials(it).(OutField),2);
    % data_trials(it).([Outxfld OutField]) = linspace(timein(1),timein(end),T);
    data_trials(it).([Outxfld OutField]) = tEst';
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
