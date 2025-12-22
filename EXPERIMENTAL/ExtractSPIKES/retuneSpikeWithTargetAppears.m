function S=retuneSpikeWithTargetAppears(S,Trials)
nAllTrials=length(Trials);
cells=size(S.AllNeuronIdentity,1);

for iT=1:nAllTrials
    for iC=1:cells
        if ~isempty(S.t{iT,iC}) && S.successfull_trials(iT)
            S.t{iT,iC}= S.t{iT,iC} - S.time_target_trials(iT);
        end
    end
end