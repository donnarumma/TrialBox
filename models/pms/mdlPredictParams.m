function   par = mdlPredictParams(parSource)
% function par = mdlPredictParams(parSource)
par.exec            = true;
par.kfold           = 4;
par.numIterations   = 100;
par.InField         = 'y';
par.OutField        = 'y';
par.SuccessField    = 'success';
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end