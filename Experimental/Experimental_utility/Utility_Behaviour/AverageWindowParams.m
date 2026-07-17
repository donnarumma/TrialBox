function   par = AverageWindowParams(parSource)
% function par = AverageWindowParams(parSource)
par.exec    = true;
par.useSqrt = false;
par.average = false;
par.binWidth= 20;
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