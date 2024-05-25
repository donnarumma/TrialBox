function   par = srssdDecodeParams(parSource)
% function par = srssdDecodeParams(parSource)
par.exec            = true;
par.InField         = 'y';
par.OutField        = 'y';
par.nChannels       = 1;
par.nTimes          = 1;
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end