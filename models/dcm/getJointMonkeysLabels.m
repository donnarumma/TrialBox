function   [lab,TargetLab]=getJointMonkeysLabels(labvec)
% function [lab,TargetLab]=getJointMonkeysLabels(labvec)

Len     = length(labvec);
CondLabs={'ACT S, OBS K', 'OBS S, ACT K','ACT S, ACT K'};
DirLabs =cell(8,1);
for iD=1:8
    DirLabs{iD}=['DIR ' num2str(iD)];
end   
if Len==8
    % 8101  Target position 1
    % 8102  Target position 2
    % 8103  Target position 3 
    % 8104  Target position 4 
    % 8105  Target position 5
    % 8106  Target position 6 
    % 8107  Target position 7 
    % 8108  Target position 8
    TargetLab=DirLabs;
end
if Len==3
    TargetLab=CondLabs;
    % 8001  SOLO cursor 1 (monkey S)
    % 8002  SOLO cursor 2 (monkey K)
    % 8004  JOINT ACTION condition 
end
if Len==3*8
    TargetLab=cell(3*8,1);
    for iC=1:3
        for iD=1:8
            TargetLab{iD + 8*(iC-1)}=[DirLabs{iD} ,', ' CondLabs{iC}];
        end
    end

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
end
lab = TargetLab{labvec};
