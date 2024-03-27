function   [data_trials,out] = SmoothWindow(data_trials, par)
% function [data_trials,out] = SmoothWindow(data_trials, par)
% 
execinfo=par.exec;
if ~isempty(execinfo); t=tic; end

InField =par.InField;
OutField=par.OutField;
binWidth=par.binWidth;

Ntrials =length(data_trials);
for iTrial = 1:Ntrials
    timeField   =data_trials(iTrial).(['time' InField]);
    tInit       =min(timeField);
    tEnd        =max(timeField);
    dt          =mean(diff(timeField));
    ti          =linspace(tInit,tEnd,length(timeField));
    nTimesOut   =length(ti)-binWidth+1;
    timeout     =(ti(1)+(binWidth+1)*dt/2):dt:(ti(end)-(binWidth-3)*dt/2);
    inData      =data_trials(iTrial).(InField);
    nChannels   =size(inData,1);
    dataout     =NaN(nChannels,nTimesOut);
    for iTimeOut=1:nTimesOut
        dataout(:,iTimeOut)=mean(inData(:,iTimeOut:iTimeOut+binWidth-1),2);
    end
    
    data_trials(iTrial).(OutField)          = dataout;
    data_trials(iTrial).(['time' OutField]) = timeout;

end
    

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end
