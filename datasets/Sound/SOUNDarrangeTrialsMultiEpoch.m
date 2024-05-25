function   [dataTrials,out]= SOUNDarrangeTrialsMultiEpoch(data_raw,par)
%function  [dataTrials,out]= SOUNDarrangeTrialsMultiEpoch(data_raw,par)
execinfo    =par.exec;
if ~isempty(execinfo); t=tic; end

InField       = par.InField;
OutField      = par.OutField;

fs          = data_raw.fsample;
N_Trials    = length(data_raw.eeg_sound);
it1         = par.it_start*fs;

trialNames  ={'Eng Coher','Ita Coher','Eng Incoher','Ita Incoher','Scrum Inchoer'};
% dataTrials     =repmat(struct('trialId',[],fname,[],'time',[],'trialType',[],'T',[]),N_Trials,1);
% dataTrials     =repmat(struct('trialId',[],fname,[],['time' fname],[],'trialType',[]),N_Trials,1);
dataTrials = struct();
for iTr = 1:N_Trials
    it2 = size(data_raw.(InField){iTr},2);
    dataTrials(iTr,1).trialId               = iTr;
    dataTrials(iTr,1).(['time' InField])    = linspace(0,size(data_raw.(InField){iTr}(1:par.ch_number,it1:it2),2)/fs,size(data_raw.(InField){iTr}(1:par.ch_number,it1:it2),2));
    dataTrials(iTr,1).events.Sound_start    = 0;
    dataTrials(iTr,1).events.Sound_end      = it2;
    dataTrials(iTr,1).(OutField)            = data_raw.(InField){iTr}(1:par.ch_number,it1:it2);
    dataTrials(iTr,1).trialType             = data_raw.y(iTr);
    dataTrials(iTr,1).trialName             = trialNames{data_raw.y(iTr)};
    dataTrials(iTr,1).train                 = false; %false
    dataTrials(iTr,1).test                  = false;    
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end