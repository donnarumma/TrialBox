function   par = plot_MontageParams(parSource)
% function par = plot_MontageParams(parSource)
par.exec        = true;
par.hfig        = [];
par.col         = [0.5,0.5,0.7]; % col borders
par.hmShow      = 49;

try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end