function   par = removeInactiveParams(parSource)
% function par = removeInactiveParams(parSource)
par.exec    = true;
par.InField = 'y';
par.OutField= 'y';
    try
        fnames=fieldnames(parSource);
        for in=1:length(fnames)
            par.(fnames{in})=parSource.(fnames{in});
        end
    catch
    end
end