function   [dataTrials,out]= AUDIOBOOKarrangeTrials(data_raw,par)
%function  [dataTrials,out]= AUDIOBOOKarrangeTrials(data_raw,par)
execinfo    =par.exec;
if ~isempty(execinfo); t=tic; end

InField       = par.InField;
OutField      = par.OutField;

fs          = par.fs;
N_Trials    = length(data_raw);

it1         = par.it_start*fs;
it2         = it1 + par.it_end * fs-1;

trialNames  ={'English','Italian'};
% dataTrials     =repmat(struct('trialId',[],fname,[],'time',[],'trialType',[],'T',[]),N_Trials,1);
% dataTrials     =repmat(struct('trialId',[],fname,[],['time' fname],[],'trialType',[]),N_Trials,1);
dataTrials = struct();
for iTr = 1:N_Trials
    dataTrials(iTr,1).trialId                   = iTr;
    dataTrials(iTr,1).(['time' InField])        = data_raw(iTr).time;
    dataTrials(iTr,1).events.Audiobook_start    = 0;
    dataTrials(iTr,1).events.Audiobook_end      = it2;
    dataTrials(iTr,1).(OutField)                = data_raw(iTr).(InField)(:,it1:it2);
    dataTrials(iTr,1).trialType                 = data_raw(iTr).trialinfo;
    dataTrials(iTr,1).trialName                 = trialNames{data_raw(iTr).trialinfo};
    dataTrials(iTr,1).train                     = false; %false
    dataTrials(iTr,1).test                      = false;    
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end