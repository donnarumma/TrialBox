function   par = plot_EachDimVsTimeParams(parSource)
% function par = plot_EachDimVsTimeParams(parSource)
par.exec        = true;
par.InField     = 'y';
par.P           = 95;      % Confidence interval 68 (1 std)- 95 (2 std)- 99.7 (3 std)
par.N           = 10000;   % number of bootstrap samples
% par.opt=[1,0,0];  % BOOTSTRAP
% par.opt=[0,1,0];  % BOOTSTRAP, equivalent to previous option - slower and demonstrative
% par.opt=[0,0,1];  % mean and variance
par.opt         = [0,0,1]; % [1,0,0] perform bootstrap vs [0,0,1] without bootstrap
par.SE          = 0; % 0 standard deviation, 1 standard error, bootstrap standard error 
par.hfig        = [];
par.cmaps       = [];
par.cmapslight  = [];
par.grayzone    = [];
par.legplot     = 2;
par.novariance  = false;
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