function     [data_trials, out] = meanData(data_trials,par)
% function   [data_trials, out] = meanData(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
OutField                = par.OutField;
InField                 = par.InField;
try
    trialTypes          = par.trialTypes;
catch
    trialTypes          = [data_trials.trialType];
end
trialNames              = {data_trials.trialName};
[classes,ord]           = unique(trialTypes);
trialNames              = trialNames(ord);
nClasses                = length(classes);
ifclass                 = par.ifclass;
cl_data_trials          = struct;
for icl=1:nClasses
    cl                  = classes(icl);
    icl_idx             = trialTypes==cl;
    dat                 = cat(3,data_trials(icl_idx).(InField));
    [nChannels,nTimes,~]= size(dat);
    mtl                 = nan(nChannels,nTimes);
    stl                 = nan(nChannels,nTimes);
    for kvar = 1:nChannels
        if nTimes==1
            data_kvar       = squeeze(dat(kvar,:,:));
        else
            data_kvar       = squeeze(dat(kvar,:,:))';
        end
        try
            [mtl(kvar,:),stl(kvar,:)]               = CI_compute(data_kvar,par);
        catch
            for iTime=1:nTimes
                [mtl(kvar,iTime),stl(kvar,iTime)]   = CI_compute(data_kvar(:,iTime),par);
            end
        end
    end
    if par.SE==1
        stl=stl/sqrt(size(data_kvar,1));
    elseif par.SE==2 && (par.opt(1) || par.opt(2))
        stl=stl/sqrt(par.N);
    end
    cl_data_trials(icl).(OutField)          = mtl;
    cl_data_trials(icl).([OutField 'sigma'])= stl;
    cl_data_trials(icl).([OutField '3d'])   = dat; % nChannels x nTimes x num Elements in class
    cl_data_trials(icl).([OutField 'nS'])   = size(dat,3);
    cl_data_trials(icl).trialName           = trialNames{icl};
    cl_data_trials(icl).trialType           = cl;
    cl_data_trials(icl).trialId             = icl;
    cl_data_trials(icl).train               = true;
    cl_data_trials(icl).valid               = false;
    cl_data_trials(icl).test                = false;
    cl_data_trials(icl).(['time' OutField]) = data_trials(1).(['time' par.InField]);
end
if ifclass
    out.data_trials =data_trials;
    data_trials     =cl_data_trials;
else
    out.data_trials =cl_data_trials;
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end