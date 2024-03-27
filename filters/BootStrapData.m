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
            cl_data_trials(icls(iSample)).(OutField)(iChannel,:) = channelBootsVals(iSample,:);
            cl_data_trials(icls(iSample)).trialName           = trialNames{icl};
            cl_data_trials(icls(iSample)).trialType           = cl;
            cl_data_trials(icls(iSample)).trialId             = icl;
            cl_data_trials(icls(iSample)).train               = true;
            cl_data_trials(icls(iSample)).valid               = false;
            cl_data_trials(icls(iSample)).test                = false;
            cl_data_trials(icls(iSample)).(['time' OutField]) = data_trials(1).(['time' par.InField]);
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