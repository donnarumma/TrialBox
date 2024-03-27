function   par = plot_scatterParams(parSource)
% function par = plot_scatterParams(parSource)
par.exec        = true;
par.hfig        = [];
par.wd          = [];
par.cmaps       = [];
par.evaltime    = [];
par.xfld        = 'time';
par.markerSize  = 10;
par.marker      = 'o';
par.filled      = true;     
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end