function   par = cspModelParams(parSource)
% function par = cspModelParams(parSource)
par.exec            = true;
par.m               = 2;
par.InField         = 'y';
par.OutField        = 'y';
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end