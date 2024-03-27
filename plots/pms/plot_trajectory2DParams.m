function   par = plot_trajectory2DParams(parSource)
% function par = plot_trajectory2DParams(parSource)
par.exec        = true;
par.hfig        = [];
par.wd          = [];
par.cmaps       = [];
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end