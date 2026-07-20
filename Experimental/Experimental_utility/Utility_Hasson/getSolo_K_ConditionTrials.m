function S=getSolo_K_ConditionTrials(S,Trials)

VARIABLES_monkey_joint_action;
nf='solo_K_trials';
tcond=SOLO_K_CONDITION;

S.(nf)=zeros(length(Trials),1);
for i=1:length(Trials)
    if ismember(tcond,Trials(i).Event(:,1));
        S.(nf)(i)=1;
    end
end