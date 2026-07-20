function S=retuneSpikeWithMovementOnsetS(S,Trials)

nAllTrials=length(Trials);
cells=size(S.AllNeuronIdentity,1);


for iT=1:nAllTrials
    if S.successfull_trials(iT)
        if S.joint_trials(iT)
%             t0=(S.time_S_movement_onset(iT)+S.time_K_movement_onset(iT))/2;
            t0=S.time_S_movement_onset(iT);
        elseif S.solo_S_trials(iT)
            t0=S.time_S_movement_onset(iT);
        elseif S.solo_K_trials(iT)
            t0=S.time_K_movement_onset(iT);
        end
        for iC=1:cells
            if ~isempty(S.t{iT,iC}) 
                S.t{iT,iC}= S.t{iT,iC} - t0;
            end
        end
    end
end