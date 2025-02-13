function   [dataTrials,out]= AUDIOBOOKarrangeTrialsMultiEpoch(data_raw,par)
%function  [dataTrials,out]= AUDIOBOOKarrangeTrialsMultiEpoch(data_raw,par)
execinfo    =par.exec;
if ~isempty(execinfo); t=tic; end

InField       = par.InField;
OutField      = par.OutField;

N_Trials    = length(data_raw);

trialNames  ={'English','Italian'};

dataTrials = struct();
for iTr = 1:N_Trials
    dataTrials(iTr,1).trialId                   = iTr;
    dataTrials(iTr,1).(['time' InField])        = data_raw(iTr).time;
    dataTrials(iTr,1).events.Audiobook_start    = data_raw(iTr).time(1);
    dataTrials(iTr,1).events.Audiobook_end      = data_raw(iTr).time(end);
    dataTrials(iTr,1).(OutField)                = data_raw(iTr).(InField);
    dataTrials(iTr,1).trialType                 = data_raw(iTr).trialinfo;
    dataTrials(iTr,1).trialName                 = trialNames{data_raw(iTr).trialinfo};
    dataTrials(iTr,1).train                     = false; %false
    dataTrials(iTr,1).test                      = false;    
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end