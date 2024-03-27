function   par = plot_trajectory3DParams(parSource)
% function par = plot_EachDimVsTimeParams(parSource)
par.exec        = true;
par.InField     = 'y';
par.hfig        = [];
par.cmaps       = [];
par.cmapslight  = [];
par.grayzone    = [];
par.keep        = 1;
par.ylabel      = '$y$';
par.nCols       = 1;
par.YLIM        = [-inf,inf];
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end