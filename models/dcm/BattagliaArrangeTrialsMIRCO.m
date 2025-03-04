function   [data_trials, session_data, Trials]= BattagliaArrangeTrialsMIRCO(par)
% function [data_trials, session_data, Trials]= BattagliaArrangeTrials(par)
% 
%% Extract a data session
selS            = par.selS;
selK            = par.selK;
dmode           = par.dmode;
session_name    = par.session_name;
session_data    = lfpSessionExtractMIRCO(session_name,dmode,selS,selK);
InField         = 'LFP';
TField          = 'tLFP';
Trials          = session_data.Trials;
nTrials         = length(Trials);
data_trials=struct;
% [~,Labels] = getJointMonkeysLabels(1:24);
for iTrial=1:nTrials
    %data_trials(iTrial).CSD         =permute(session_data.CSD{iTrial},[2,3,1]);
    %data_trials(iTrial).(['H' 'CSD']) =session_data.Hz;
    % data_trials(iTrial).(TField)    =session_data.Trials(iTrial).RawTime(2:end);
    data_trials(iTrial).(TField)    =session_data.Trials(iTrial).RawTime;
    data_trials(iTrial).(InField)   =session_data.LFP{iTrial}';
    data_trials(iTrial).JXYEXY      =session_data.Trials(iTrial).JXYEXY;
    data_trials(iTrial).timeJXYEXY  =session_data.Trials(iTrial).timeJXYEXY;
    data_trials(iTrial).dt          =session_data.dt;
    data_trials(iTrial).selIndex    =session_data.selIndex;
    data_trials(iTrial).iLFP_K      =session_data.iLFP_K;
    data_trials(iTrial).iLFP_S      =session_data.iLFP_S;
    data_trials(iTrial).ARp         =session_data.ARp;
    data_trials(iTrial).rsfactor    =session_data.rsfactor;
    data_trials(iTrial).AlignEvent  =session_data.Trials(iTrial).AlignEvent;
    data_trials(iTrial).Chamber     =session_data.Trials(iTrial).Chamber;
    data_trials(iTrial).session_name=par.session_name;
    data_trials(iTrial).trialId     =iTrial;
    data_trials(iTrial).trialType   =find(session_data.klabels(iTrial,:));
    % data_trials(iTrial).trialName   =Labels{data_trials(iTrial).trialType};
    data_trials(iTrial).klabels     =session_data.klabels(iTrial,:);
    data_trials(iTrial).RT_S        =session_data.Trials(iTrial).RT_S;
    data_trials(iTrial).RT_K        =session_data.Trials(iTrial).RT_K;
    data_trials(iTrial).ET          =session_data.Trials(iTrial).ET;
    data_trials(iTrial).timeET      =session_data.Trials(iTrial).timeET;
end

