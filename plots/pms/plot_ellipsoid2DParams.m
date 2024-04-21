function   par = plot_ellipsoid2DParams(parSource)
% function par = plot_ellipsoid2DParams(parSource)
par             = CI_computeParams;
par.exec        = true;
par.hfig        = [];
par.cmaps       = [];
par.InField     = 'xorth';
par.nCols       = 2;
par.nRows       = 2;

try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end