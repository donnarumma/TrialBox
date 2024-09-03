function [x0s, labs, labnames] = getPointsPCA(seq,selected,indT,dimsToPlot,xspec)
Lseq=length(seq);
% x0s=(lseq);
% labs=[];
try
    xspec;
catch
    xspec='xorth';
end
labs=[];
x0s =[];
labnames=cell(0,0);
count=0;
for n=1:Lseq
    if ismember(seq(n).trialType,selected)
        count           = count+1;
        labs(:,count)   = seq(n).trialType;
        labnames{count} = seq(n).trialName;
        x0s(count,:)    = seq(n).(xspec)(dimsToPlot,indT);
    end
end
return
%% old version: warning no nan is possible
Lseq=length(seq);
x0s=[];
labs=[];
xspec='xorth';
for n=1:Lseq
    if ismember(seq(n).trialType,selected)
        labs=[labs,seq(n).trialType];
        x0s = [x0s; seq(n).(xspec)(dimsToPlot,indT)];
    end
end
return   
[~,b]=find(~isnan(x0s));
b=unique(b);
    
dL=Lseq-length(b);
if dL>0
    fprintf('Warning: removed %g columns\n',dL);
end
    
x0s=x0s(:,b);
labs=labs(b);
    
x0s=x0s';
