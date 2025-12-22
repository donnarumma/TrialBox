function S=getTimeTargetTrials(S,Trials)

VARIABLES_monkey_joint_action;
nf='time_target_trials';
tcond=TARGET_APPEARS;

S.(nf)=nan(length(Trials),1);
for i=1:length(Trials)

    if ismember(SUCCESS,Trials(i).Event(:,1));
        v=Trials(i).Event(:,1)==tcond;
        S.(nf)(i)=Trials(i).Event(v,2);
    end
end