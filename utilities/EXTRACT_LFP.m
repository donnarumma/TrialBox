function LFP_H=EXTRACT_LFP(Raw_H,Raw_E,params)
% function LFP_H=EXTRACT_LFP(Raw_H,Raw_E)
%
%     % resampling parameters
%     params.RS_Factor        =20;  % resampling factor   
%
%     % apply filter params
%     params.filterOrder      =100;
%
%     % NOTCH params
%     params.NOTCHfilter      =0;   % apply NOTCH filter
%     params.NOTCHfilterLow   =45;  % low   NOTCH freq
%     params.NOTCHfilterHigh  =55;  % high  NOTCH freq
%
%     params.LPfilter         =0;   % apply   LP  filter
%     params.LPfiilterLow     =1;   % low     LP  freq
%     params.LPfilterHigh     =100; % high    LP  freq
%
%     params.SubtractMean     =0;   % subtract mean values
%
%     params.SubtractDCoffset =1;   % subtract DC OFFSET
%
%     params.removeEYE        =0;   % remove EYE artifacts
%
%     params.removeRAW        =1;   % remove field RAW after extraction
    
default_params=getDefaultLFPparams;
try
    params=recopyFields(params,default_params);
catch
    params=default_params;
end
    ti = tic;
    LFP_H   = Raw_H;
    LFP_H   = Odysseas_Resample_Trials_params(LFP_H,params);
    LFP_H   = applyFilterChain(LFP_H,params);
    
    if params.removeEYE
        LFP_E   = Raw_E;
        LFP_E   = Odysseas_Resample_Trials_params(LFP_E,params);
        LFP_E   = applyFilterChain(LFP_E,params);
        LFP_H   = Odysseas_Get_Eye_Artifact_debug(LFP_H,LFP_E);
%         LFP_H   = Odysseas_Get_Eye_Artifact(LFP_H,LFP_E);
    end
    if params.removeRAW
        LFP_H   = rmfield(LFP_H,'Raw');
    end
    fprintf('Elapsed time: %g s\n',toc(ti));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function Trials=applyFilterChain(Trials,params)
%% function Trials=applyFilterChain(Trials,params)

    f_Ny    =Trials(1).FS/2;
    EndCh   =size(Trials(1).Raw,1);
    nTrials =length(Trials);
    order   =params.filterOrder;
    
    try 
        Session = Trials(1).Session;
        fprintf('Filtering LFP data from %s\n',Session);
    catch
    end

    if params.NOTCHfilter
        % NOTCH filter set up
        n1=params.NOTCHfilterLow/f_Ny;
        n2=params.NOTCHfilterHigh/f_Ny;
        NT_b = fir1(order,[n1 n2],'stop');
    end
    if params.LPfilter
        % LP filter set up
        q1=1/f_Ny;
        q2=100/f_Ny;
        LP_bq = fir1(order,[q1 q2]);              
    end

    for iT=1:nTrials
       % loop on all channels
       for channel=1:EndCh
            Trial                       = double(Trials(iT).Raw(channel,:));
            if params.NOTCHfilter
                % apply NOTCH filter
                Trial                       = filtfilt(NT_b,1,Trial);
            end
            if params.LPfilter
                % apply LP filter
                Trial                       = filtfilt(LP_bq,1,Trial);  
            end
            if params.SubtractMean
                % subtract mean value
                Trial                       = demean(Trial);
            end
            if params.SubtractDCoffset
                % subtract DC OFFSET
                Trial                       = detrend(Trial);
            end        
            Trials(iT).LFP(channel,:)       = Trial;
       end
    end
end