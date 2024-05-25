function   par = cspEncodeParams(parSource)
% function par = cspEncodeParams(parSource)
par.exec        = true;
par.W           = [];
par.InField     = 'y';
par.OutField    = 'y';
par.xfld        = 'time';
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end