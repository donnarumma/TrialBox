function   par = dataSplitParams(parSource)
% function par = dataSplitParams(parSource)
par.exec            = true;
par.TrainPercentage = 70;
par.TrainPercentageV= 0;
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