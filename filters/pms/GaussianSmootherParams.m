function   par = GaussianSmootherParams(parSource)
% function par = GaussianSmootherParams(parSource)
par.exec    = 1;
par.causal  = false;    % logical indicating whether temporal smoothing shouldw include only past data (true) or all data (false)
par.kernSD  = 30;       % standard deviation of Gaussian kernel, in msec
par.stepSize= 20;       % time between 2 consecutive datapoints in yIn, in msec
par.nSD     = 3;        % threshold of nSD standard deviations out
par.InField = 'y';
par.OutField= 'y';
    try
        fnames=fieldnames(parSource);
        for in=1:length(fnames)
            par.(fnames{in})=parSource.(fnames{in});
        end
    catch
    end
end