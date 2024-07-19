function   [dataTrials,out]= SOUNDarrangeTrialsMultiEpoch(data_raw,par)
%function  [dataTrials,out]= SOUNDarrangeTrialsMultiEpoch(data_raw,par)
execinfo    =par.exec;
if ~isempty(execinfo); t=tic; end

InField       = par.InField;
OutField      = par.OutField;

N_Trials    = length(data_raw.eeg);

trialNames  ={'Eng Coher','Ita Coher','Eng Incoher','Ita Incoher','Scrum Inchoer'};
% dataTrials     =repmat(struct('trialId',[],fname,[],'time',[],'trialType',[],'T',[]),N_Trials,1);
% dataTrials     =repmat(struct('trialId',[],fname,[],['time' fname],[],'trialType',[]),N_Trials,1);
Tab = data_raw.Tab_final;
N = size(data_raw.(InField){1},2);
y_approx = 9;
delta =  y_approx / (N - 1500);
x = -1499 * delta;
time = linspace(x, y_approx, N);
dataTrials = struct();
for iTr = 1:N_Trials
    dataTrials(iTr,1).trialId               = iTr;
    dataTrials(iTr,1).(['time' InField])    = time;
    dataTrials(iTr,1).events.Text_start     = time(Tab(iTr).StartText);
    dataTrials(iTr,1).events.Text_end       = time(Tab(iTr).EndText);
    dataTrials(iTr,1).events.Sound_start    = time(Tab(iTr).StartSound);
    dataTrials(iTr,1).events.Sound_end      = time(Tab(iTr).EndSound);
    dataTrials(iTr,1).(OutField)            = data_raw.(InField){iTr};
    dataTrials(iTr,1).trialType             = data_raw.trialinfo(iTr);
    dataTrials(iTr,1).trialName             = trialNames{data_raw.trialinfo(iTr)};
    dataTrials(iTr,1).train                 = false; %false
    dataTrials(iTr,1).test                  = false;    
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end