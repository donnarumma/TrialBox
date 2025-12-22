function S=prepareSpikeAndTrials(S,Trials,interval,selMonkey)

try
    interval;
catch
    interval=[-100,300];
    interval=[-500,1000];
end

try
    selMonkey;
catch
    selMonkey='all';
end
if strcmp(selMonkey,'all')
    selMonkey =  '1';     % all
elseif strcmp(selMonkey,'S')
    selMonkey =  'iC<6'; % -> S
elseif strcmp(selMonkey,'K')
%    selMonkey =  'iC>4'; % -> K
    selMonkey =  'iC>5'; % -> K
end

%start_time=-100;end_time=300;
start_time=min(interval);
end_time  =max(interval);

nAllTrials=length(Trials);
cells=size(S.AllNeuronIdentity,1);
bin=1/1000;
time=start_time/1000:bin:(end_time/1000+bin);
joint_trials =S.joint_trials;
solo_S_trials=S.solo_S_trials*2;
solo_K_trials=S.solo_K_trials*3;
trialType2=joint_trials+solo_S_trials+solo_K_trials;

trialType=zeros(nAllTrials,1);
for i=1:8
    trialType=trialType+S.(['T' num2str(i) '_trials']) * i;
end


good_iCs=[];
for iC=1:cells
    if eval(selMonkey)
        good_iCs=[good_iCs iC];
    end
end

good_cells=length(good_iCs);
count=0;
for iT=1:nAllTrials
    if S.successfull_trials(iT)
        sp=zeros(good_cells,length(time));
        count=count+1;
        for inC=1:length(good_iCs)
            iC=good_iCs(inC);
            if ~isempty(S.t{iT,iC})
                sp(inC,:)=hist(S.t{iT,iC},time);
            end
        end
        dat(count).trialId   =iT;
        dat(count).spikes    =sp(:,2:end-1);
        dat(count).trialType =trialType(iT);
        dat(count).trialType2=trialType2(iT);
    end
end
S.dat=dat;