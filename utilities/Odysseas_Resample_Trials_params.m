function Trials=Odysseas_Resample_Trials_params(Trials,params)
% function Trials=Odysseas_Resample_Trials_params(Trials,params)    
    if round(Trials(1).FS)==round(2.441406250000000e+04)
        RS_Factor=params.RS_Factor;
    %     RS_Factor=20;
    else
        RS_Factor=10;
        fprintf('Warning: Resampling=%g',RS_Factor); pause;
    end

    New_FS=Trials(1).FS/RS_Factor;
    EndCh=size(Trials(1).Raw,1);
      
    %Remap
    for i=1:size(Trials,2)
        Old_Raw     =Trials(i).Raw;
        Old_InterRaw=Trials(i).InterRaw;
        
        %Preallocate       
        New_End     =round(size(Trials(i).Raw,2)/RS_Factor);
        Trials(i).Raw(:,New_End+1:end)=[];           
        New_End     =round(size(Trials(i).InterRaw,2)/RS_Factor);
        Trials(i).InterRaw(:,New_End+1:end)=[];   
           
        for channel=1:EndCh
               
            %Signal
            x=Old_Raw(channel,:);
            NANs=RS_Factor*round(size(x,2)/RS_Factor)-size(x,2);
            if NANs>0 
               x(end+1:size(x,2)+NANs)=zeros(NANs,1);
            elseif NANs<0
               x(end+1+NANs:end)=[];
            else
               %Nothing
            end        
            x=reshape(x,RS_Factor,[]);       
            x=mean(x);
            Trials(i).Raw(channel,:)=x;
            
            %InterTrial
            y=Old_InterRaw(channel,:);
            NANs=RS_Factor*round(size(y,2)/RS_Factor)-size(y,2);
            if NANs>0 
               y(end+1:size(y,2)+NANs)=zeros(NANs,1);
            elseif NANs<0
               y(end+1+NANs:end)=[];
            else
               %Nothing
            end        
            y=reshape(y,RS_Factor,[]);
            y=mean(y);       
            Trials(i).InterRaw(channel,:)=y;
           
       end
       
       %TimeStamp
       Trials(i).TimeStamp=[Trials(i).TimeStamp(1),Trials(i).TimeStamp(1)+size(x,2)/New_FS];
       
       %FS
       Trials(i).FS=New_FS;
      
       %RawTime
        T=Trials(i).RawTime;
        NANs=RS_Factor*round(size(T,2)/RS_Factor)-size(T,2);
        if NANs>0 
            T(end+1:size(T,2)+NANs)=zeros(NANs,1);
        elseif NANs<0
            T(end+1+NANs:end)=[];
        else
            %Nothing
        end        
        T=reshape(T,RS_Factor,[]);
        Trials(i).RawTime=T(1,:);       
    
    end
  
    clearvars -except Trials

end