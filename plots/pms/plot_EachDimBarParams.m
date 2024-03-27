function   par = plot_EachDimBarParams(parSource)
% function par = plot_EachDimBarParams(parSource)
par.exec        = true;
par.hfig        = [];
par.grayzone    = [];
par.novariance  = false;
par.addbar      = false;
par.cmaps       = [];
par.legplot     = 2;
par.nCols       = 1;
par.ylabel      = '$y$';
par.evaltime    = [];%1;
par.chanceline  = false;
par.cmaps       = [];
par.cmapslight  = [];
par.grayzone    = [];
par.legplot     = 2;
par.novariance  = false;
par.keep        = 1;
par.YLIM        = [-inf,inf];
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end