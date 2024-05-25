function   par = fitNBPWParams(parSource)
% function par = fitNBPWParams(parSource)
par.exec            = true;
    try
        fnames=fieldnames(parSource);
        for in=1:length(fnames)
            par.(fnames{in})=parSource.(fnames{in});
        end
    catch
    end
end