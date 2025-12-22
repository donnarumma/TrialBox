function S=getTConditionTrials(S,Trials,ind)
try
    ind;
catch
    ind=1;
end
VARIABLES_monkey_joint_action;
nf    =['T' num2str(ind) '_trials'];
tcond =eval(['T' num2str(ind)]);

S.(nf)=zeros(length(Trials),1);
for i=1:length(Trials)
    if ismember(tcond,Trials(i).Event(:,1));
        S.(nf)(i)=1;
    end
end