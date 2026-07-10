function Trials = Donnarumma_Align_Trials(Trials, AlignEvent,Start,Stop)
%Re-alignes trials to the selected event %Trials=trial structure, Al
% es. Trials= Realign(Trials, 8001,-0.7,0.7)
goodtrials=true(length(Trials),1);
for i=1:size(Trials,2)  
    TimeLock=Trials(i).Event (Trials(i).Event (:,1)==AlignEvent,2); %Get relative time of event presentation    
    if length(TimeLock)>1;TimeLock=TimeLock(end);end %Multiple 28 events
    if ~isempty(TimeLock)
        [StartP,StartT] =Odysseas_Get_Nearest_Value(Start+TimeLock ,Trials(i).RawTime);%Get start centered on the event
        [StopP,StopT] =Odysseas_Get_Nearest_Value(Stop+TimeLock,Trials(i).RawTime);  %Get stop centered on the event      
        

        %Select the field
        if isfield(Trials, 'Raw')
            Trials(i).Raw=Trials(i).Raw(:,StartP:StopP); %Get the newly centered matrix
        elseif isfield(Trials, 'LFP')
            Trials(i).LFP=Trials(i).LFP(:,StartP:StopP); %Get the newly centered matrix
        elseif isfield(Trials, 'MUA')
            Trials(i).MUA=Trials(i).MUA(:,StartP:StopP); %Get the newly centered matrix
        elseif isfield(Trials, 'TFC')
            Trials(i).TFC=Trials(i).TFC(:,StartP:StopP,:); %Get the newly centered matrix 
        elseif isfield(Trials, 'SUA')
            Trials(i).SUA=Trials(i).SUA(StartP:StopP,:); %Get the newly centered matrix
        elseif isfield(Trials, 'TFXC')
            for channel_1=1:size(Trials(i).TFXC,1)
                for channel_2=1:size(Trials(i).TFXC,2)
                    Trials(i).TFXC(channel_1,channel_2).WC=Trials(i).TFXC(channel_1,channel_2).WC(:,StartP:StopP);
                    isfield(Trials(1).TFXC,'S1');Trials(i).TFXC(channel_1,channel_2).S1=Trials(i).TFXC(channel_1,channel_2).S1(:,StartP:StopP);
                    isfield(Trials(1).TFXC,'S2');Trials(i).TFXC(channel_1,channel_2).S2=Trials(i).TFXC(channel_1,channel_2).S2(:,StartP:StopP);
                end
            end
         elseif  isfield(Trials,'TFWC')
            for channel_1=1:size(Trials(i).TFWC,1)
                for channel_2=1:size(Trials(i).TFWC,2)
                    Trials(i).TFWC(channel_1,channel_2).WC=Trials(i).TFWC(channel_1,channel_2).WC(:,StartP:StopP);
                end
            end
        end           
        
        Trials(i).RawTime=Trials(i).RawTime(StartP:StopP);%Redefine time
        Trials(i).AlignEvent=AlignEvent; % Just in case.. 
        
        if isfield(Trials, 'JXYEXY')
            NewTime=linspace(Trials(i).TimeStamp(1),Trials(i).TimeStamp(2),size(Trials(i).JXYEXY,2));
            [StartP,~] =Odysseas_Get_Nearest_Value(Start+TimeLock ,NewTime);%Get start centered on the event
            [StopP,~] =Odysseas_Get_Nearest_Value(Stop+TimeLock,NewTime);  %Get stop centered on the event      
            Trials(i).JXYEXY=Trials(i).JXYEXY(:,StartP:StopP); %Get the newly centered matrix
            Trials(i).timeJXYEXY = NewTime(StartP:StopP);
        end
               
        
        Trials(i).TimeStamp= [StartT,StopT]; %Redefine timestamp
    else
        fprintf('Warning: Event not present in Trial %g\n',i);
        goodtrials(i)=false;
    end 
end
Trials=Trials(goodtrials);
end

