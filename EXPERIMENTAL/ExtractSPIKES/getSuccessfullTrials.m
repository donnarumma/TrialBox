function S=getSuccessfullTrials(S,Trials);

VARIABLES_monkey_joint_action;
nf='successfull_trials';
tcond=SUCCESS;
S.(nf)=zeros(length(Trials),1);
for i=1:length(Trials)
    if ismember(tcond,Trials(i).Event(:,1));
        S.(nf)(i)=1;
    end
end