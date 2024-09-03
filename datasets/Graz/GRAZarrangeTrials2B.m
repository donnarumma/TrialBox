function   [dataTrials,out]= GRAZarrangeTrials2B(data_raw,par)
%function  [dataTrials,out]= GRAZarrangeTrials2B(data_raw,par)
execinfo    =par.exec;
if ~isempty(execinfo); t=tic; end
Infield       =par.InField;
N_Trials    =length(data_raw);
it1         =par.it1;
it2         =par.it2;
fs          =par.fsample;
it_zero     =abs(it1)*fs;
n_samp      =(it2-it1)*fs;
trialNames  ={'Left Hand','Right Hand'};
try label = [data_raw.label];
catch
    label = [data_raw.truelab];
end
dataTrials = struct();
for iTr = 1:N_Trials
    dataTrials(iTr,1).trialId                       = iTr;
    dataTrials(iTr,1).(['time' Infield])            = linspace(it1,it2,n_samp);
    dataTrials(iTr,1).(['time' Infield])(it_zero)   = 0;
    dataTrials(iTr,1).events.FixationCross          = -3;
    dataTrials(iTr,1).events.Cue                    = 0;
    dataTrials(iTr,1).events.MotorImag_start        = 1;
    dataTrials(iTr,1).events.MotorImag_end          = 4;
    dataTrials(iTr,1).events.End                    = size(data_raw(iTr).(Infield));
    dataTrials(iTr,1).(Infield)                     = data_raw(iTr).(Infield);
    dataTrials(iTr,1).trialType                     = label(iTr);
    dataTrials(iTr,1).trialName                     = trialNames{label(iTr)};
    dataTrials(iTr,1).train                         = true;
    dataTrials(iTr,1).test                          = false;    
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end