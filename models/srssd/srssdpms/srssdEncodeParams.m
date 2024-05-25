function   par = srssdEncodeParams(parSource)
% function par = srssdModelParams(parSource)
par.exec            = true;
par.InField         = 'y';
par.OutField        = 'y';
par.Wsrssd          = [];           % nAtoms x nChannels*nTimes
par.xfld            = 'time';

try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end