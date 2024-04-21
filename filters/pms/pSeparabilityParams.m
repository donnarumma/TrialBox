function   par = pSeparabilityParams(parSource)
% function par = pSeparabilityParams(parSource)
par.exec        = true;
par.InField     = 'y';
par.OutField    = 'comparisons';
par.xfld        = 'time';
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end