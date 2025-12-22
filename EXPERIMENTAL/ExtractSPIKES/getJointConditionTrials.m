function S=getJointConditionTrials(S,Trials)

VARIABLES_monkey_joint_action;
nf='joint_trials';
tcond=JOINT_CONDITION;

S.(nf)=zeros(length(Trials),1);
for i=1:length(Trials)
    if ismember(tcond,Trials(i).Event(:,1));
        S.(nf)(i)=1;
    end
end