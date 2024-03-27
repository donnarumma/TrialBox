function   par = plot_AccuracyBarsParams(parSource)
% function par = plot_AccuracyBarsParams(parSource)
par.exec        = true;
par.hfig        = [];
par.cmaps       = [];
par.evaltime    = [];%0%1;
par.nCols       = 1;
par.explained   = [];
par.channelSets = [];
par.SuccessField='success';
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end