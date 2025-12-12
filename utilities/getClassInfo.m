function   [condition, direction, klabels, success] = getClassInfo(Trials)
% function [condition, direction, success] = getClassInfo(Trials)
% 8101  Target position 1
% 8102  Target position 2
% 8103  Target position 3 
% 8104  Target position 4 
% 8105  Target position 5
% 8106  Target position 6 
% 8107  Target position 7 
% 8108  Target position 8
% 3000  Trial Success - Arrive at Target (periferic) + stay in target (variable time)
% 8001  SOLO cursor 1 (monkey S)
% 8002  SOLO cursor 2 (monkey K)
% 8004  JOINT ACTION condition 
% klabels <-  condition x direction
%       1 <-       1         1  
%       2 <-       1         2  
%       3 <-       1         3  
%       4 <-       1         4  
%       5 <-       1         5  
%       6 <-       1         6  
%       7 <-       1         7  
%       8 <-       1         8  
%       9 <-       2         1  
%      10 <-       2         2  
%      11 <-       2         3  
%      12 <-       2         4  
%      13 <-       2         5  
%      14 <-       2         6  
%      15 <-       2         7  
%      16 <-       2         8  
%      17 <-       3         1  
%      18 <-       3         2  
%      19 <-       3         3  
%      20 <-       3         4  
%      21 <-       3         5  
%      22 <-       3         6  
%      23 <-       3         7  
%      24 <-       3         8  


Conds       = [8001,8002,8004];
Target      = 8101:8108;
NConds      = length(Conds);
NTargets    = length(Target);
NTrials     = length(Trials);
success     = false(NTrials,1);
direction   = false(NTrials,NTargets);
condition   = false(NTrials,NConds); 
klabels     = false(NTrials,NTargets*NConds);
for iTrial = 1:NTrials
    ev=Trials(iTrial).Event(:,1);
    success(iTrial)=sum(ev==3000)>0;
    for id=1:NTargets
        direction(iTrial,id)=sum(ev==Target(id))>0;
    end
    for ic=1:NConds
        condition(iTrial,ic)=sum(ev==Conds(ic))>0;
    end
    klabels(iTrial,:)=kron(condition(iTrial,:),direction(iTrial,:));
end
