function   par = SmoothWindowParams(parSource)
% function par = SmoothWindowParams(parSource)
par.exec    = true;
par.useSqrt = false;
par.average = false;
par.binWidth= 20;    % bins over apply moving window average
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