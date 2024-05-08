function   [dataTrials,out]= FITTSarrangeTrials(data_raw,par)
%function  [dataTrials,out]= FITTSarrangeTrials(data_raw,par)
execinfo    =par.exec;
if ~isempty(execinfo); t=tic; end
fname       =par.InField;
N_Trials    =length(data_raw.trial);
it1         =par.it_start;
it2         =par.it_end;
event       =par.event;
trialNames  ={'Small Left','Small Right','Big Left','Big Right'};
% dataTrials     =repmat(struct('trialId',[],fname,[],'time',[],'trialType',[],'T',[]),N_Trials,1);
% dataTrials     =repmat(struct('trialId',[],fname,[],['time' fname],[],'trialType',[]),N_Trials,1);
dataTrials = struct();
Tab = EventExtract(data_raw.trialinfo,event);
for iTr = 1:N_Trials
    dataTrials(iTr,1).trialId           = iTr;
    dataTrials(iTr,1).(['time' fname])  = round(data_raw.time{iTr}(it1:it2),3);
    dataTrials(iTr,1).events.Start             = Tab(iTr).Start;
    dataTrials(iTr,1).events.TargetON          = Tab(iTr).TargetON;
    dataTrials(iTr,1).events.TargetOFF         = Tab(iTr).TargetOFF;
    dataTrials(iTr,1).events.GO                = Tab(iTr).GO;
    dataTrials(iTr,1).events.MovInit           = Tab(iTr).GO + Tab(iTr).RT;
    dataTrials(iTr,1).events.Touch             = Tab(iTr).Touch;    
    dataTrials(iTr,1).(fname)           = data_raw.trial{iTr}(1:par.ch_number,it1:it2);
    dataTrials(iTr,1).trialType         = data_raw.trialinfo(iTr,11);
    dataTrials(iTr,1).trialName         = trialNames{data_raw.trialinfo(iTr,11)};
    dataTrials(iTr,1).train             = true;
    dataTrials(iTr,1).test              = false;
    % dataTrials(iTr,1).T          = length(dataTrials(iTr,1).time);

    
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end