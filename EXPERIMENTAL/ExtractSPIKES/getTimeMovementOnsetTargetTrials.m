function S=getTimeMovementOnsetTargetTrials(S,Trials)

VARIABLES_monkey_joint_action;

nf='time_S_movement_onset';
tcond=S_MOVEMENT_ONSET;
S.(nf)=nan(length(Trials),1);
for i=1:length(Trials)
    if S.joint_trials(i) || S.solo_S_trials(i)
        if ismember(SUCCESS,Trials(i).Event(:,1));
            v=Trials(i).Event(:,1)==tcond;
            S.(nf)(i)=Trials(i).Event(v,2);
        end
    end
end

nf='time_K_movement_onset';
tcond=K_MOVEMENT_ONSET;
S.(nf)=nan(length(Trials),1);
for i=1:length(Trials)
    if S.joint_trials(i) || S.solo_K_trials(i)
        if ismember(SUCCESS,Trials(i).Event(:,1));
            v=Trials(i).Event(:,1)==tcond;
            S.(nf)(i)=Trials(i).Event(v,2);
        end
    end
end