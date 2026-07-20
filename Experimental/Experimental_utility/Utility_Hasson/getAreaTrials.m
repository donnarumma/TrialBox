function S=getAreaTrials(S,Trials)

VARIABLES_monkey_joint_action;
nf='area';
tcond=SOLO_S_CONDITION;
chamb=CHAMBERS{Trials(1).chamberID(1)};
S.(nf)=chamb;
