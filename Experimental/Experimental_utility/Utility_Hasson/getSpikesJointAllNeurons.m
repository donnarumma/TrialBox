function [S, AllNeuronIdentity]=getSpikesJointAllNeurons(Trials)
%function getSpikesJOINT(data,nf)
try
    Trials;
catch
    nf='SK033H_SUA.mat';
    d=load (nf);
    Trials=d.Trials;
end

nAllTrials=length(Trials);
dt=10/1000;

NeuronIdentityK=[];
NeuronIdentityS=[];

try
    whichTrials;
catch    
    whichTrials=1:nAllTrials;
end

AllNeuronIdentity=[];
interval=nan(nAllTrials,1);

for iT=1:nAllTrials
    SUA=Trials(iT).SUA;
%    ng=~((SUA(:,3)==0) | (SUA(:,3)==9) | (SUA(:,3)==20));
    ng=((SUA(:,3)>0) & (SUA(:,3)<9));
    SUA=SUA(ng,:);
    NeuronIdentity=unique(SUA(:,2:3),'rows');
    AllNeuronIdentity=[AllNeuronIdentity;NeuronIdentity];
    interval(iT)=(Trials(iT).TimeStamp(2)-Trials(iT).TimeStamp(1));
end

AllNeuronIdentity=unique(AllNeuronIdentity,'rows');
cells=size(AllNeuronIdentity,1);
%cells=size(NeuronIdentity,1);
S.interval=interval;
S.tvec=0:dt:max(interval);
S.t=cell(nAllTrials,cells);
%S.t_abs=cell(nAllTrials,cells);

for iT=1:nAllTrials

    SUA=Trials(iT).SUA;
%    ng=~((SUA(:,3)==0) | (SUA(:,3)==9) | (SUA(:,3)==20));
    ng=((SUA(:,3)>0) & (SUA(:,3)<9));
    SUA=SUA(ng,:);    
%    S.tvec=1:dt:(Trials(iT).TimeStamp(2)-Trials(iT).TimeStamp(1));
    
%    goodTimes=size(SUA,1);
%    NeuronIdentity=unique(SUA(:,2:3),'rows');
    A=SUA(:,2:3);
%    wo y= sum(repmat(NeuronIdentity(iC,:),size(AllNeuronIdentity,1),1)==AllNeuronIdentity,2)==2;
    
    for iC=1:cells
        vv            =sum(repmat(AllNeuronIdentity(iC,:),size(A,1),1)==A,2)==2;
        if sum(vv)>0
%            S.t{iT,iC}    =SUA(vv,1)' - Trials(iT).TimeStamp(1);
            S.t{iT,iC}=SUA(vv,1)';
        else
            S.t{iT,iC}    =[];
%            S.t{iT,iC}=[];
        end
        S.label{iT,iC}=['e' num2str(AllNeuronIdentity(iC,1)) 't' num2str(AllNeuronIdentity(iC,2))];
    end
%    S.monkeyK(iT,1)=sum(NeuronIdentity(:,1)>5);
%    S.monkeyS{iT,1}=sum(NeuronIdentity(:,1)<6);
    
%    if    SUA(t,2)>5
%        Monkey='K';
%        NeuronIdentityK=[NeuronIdentityK;[SUA(t,2:3)]];
%    else
%        Monkey='S';
%        NeuronIdentityS=[NeuronIdentityS;[SUA(t,2:3)]];
%    end
%    fprintf('Electrode: %g (Monkey:%s), Tag: %g, Spike Time= %gs\n',SUA(t,1),Monkey,SUA(t,3),SUA(t,1))
end
S.AllNeuronIdentity=AllNeuronIdentity;
return
cfg.freqTheta   = 8;
cfg.ifmean      = 0;
cfg.trialsToPlot=10;
cfg.newfigure   = 1;
PlotSpikeRaster_mod(cfg,S);
hold on;
if length(cfg.trialsToPlot)<2
    nn=AllNeuronIdentity(:,1)<6;
    lnn=sum(nn)+0.5;
    xl=get(gca,'xlim');
    plot(xlim,[lnn,lnn],'r--')
    xlim(xl);
end

%dy=0.6;
%yl = [min(whichTrials)-dy, max(whichTrials)+dy];
%ylim(yl);