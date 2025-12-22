function S=getSolo_S_ConditionTrials(S,Trials)

VARIABLES_monkey_joint_action;
nf='solo_S_trials';
tcond=SOLO_S_CONDITION;

S.(nf)=zeros(length(Trials),1);
for i=1:length(Trials)
    if ismember(tcond,Trials(i).Event(:,1));
        S.(nf)(i)=1;
    end
end