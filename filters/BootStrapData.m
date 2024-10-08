function     [data_trials, out] = BootStrapData(data_trials,par)
% function   [data_trials, out] = BootStrapData(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
OutField                = par.OutField;
InField                 = par.InField;
try
    trialTypes          = par.trialType;
catch
    trialTypes          = [data_trials.trialType];
end
trialNames              = {data_trials.trialName};
[classes,ord]           = unique(trialTypes);
trialNames              = trialNames(ord);
nClasses                = length(classes);
ifclass                 = par.ifclass;
P                       = par.P;
nBootSamples            = par.N;

cl_data_trials          = struct;
alpha                   =(1-P/100);
iclend                  = 0;
for icl=1:nClasses
    cl                  = classes(icl);
    icl_idx             = trialTypes==cl;
    dat                 = cat(3,data_trials(icl_idx).(InField));
    [nChannels,nTimes,~]= size(dat);
    iclstart = iclend+1;
    iclend   = nBootSamples*(icl);
    icls = iclstart:iclend;
    if par.TimeLagged
        % add lagged trials
        tmpdat  = dat;        
        for ilag    = 1:par.TimeLagged
            adddat  = tmpdat;
            adddat(:,1:(end-ilag),:)=tmpdat(:,(1+ilag):end,:);
            dat     = cat(3,dat,adddat);
            adddat  = tmpdat;        
            adddat(:,(1+ilag):end,:)=tmpdat(:,1:(end-ilag),:);
            dat     = cat(3,dat,adddat);
        end
    end

    for iChannel = 1:nChannels      
        if nTimes==1
            data_channel       = squeeze(dat(iChannel,:,:));
        else
            data_channel       = squeeze(dat(iChannel,:,:))';
        end
        channelBootsVals       = nan(nBootSamples,nTimes);
        for iTime=1:nTimes
            % [mtl(kvar,iTime),stl(kvar,iTime)]   = CI_compute(data_channel(:,iTime),par);
            [~,channelBootsValsTime] = bootci(nBootSamples,{@mean,data_channel(:,iTime)},'alpha',alpha);
            channelBootsVals(:,iTime)=channelBootsValsTime;
        end
        
        
        for iSample=1:nBootSamples
            cl_data_trials(icls(iSample),1).(OutField)(iChannel,:) = channelBootsVals(iSample,:);
            cl_data_trials(icls(iSample),1).trialName           = trialNames{icl};
            cl_data_trials(icls(iSample),1).trialType           = cl;
            cl_data_trials(icls(iSample),1).trialId             = icl;
            cl_data_trials(icls(iSample),1).train               = true;
            cl_data_trials(icls(iSample),1).valid               = false;
            cl_data_trials(icls(iSample),1).test                = false;
            cl_data_trials(icls(iSample),1).(['time' OutField]) = data_trials(1).(['time' par.InField]);
        end
    end
end
if ifclass
    out.data_trials =data_trials;
    data_trials     =cl_data_trials;
else
    out.data_trials =cl_data_trials;
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end