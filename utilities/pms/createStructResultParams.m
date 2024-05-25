function   par = createStructResultParams(parSource)
% function par = createStructResultParams(parSource)
par.exec            = true;
% Filter
par.attenuation     = 40;
% Model
par.numIterations   = 100;
par.kfold           = 4;
par.cspModel.m      = 0;
par.InField         = 'y';
par.OutField        = 'y';
    try
        fnames=fieldnames(parSource);
        for in=1:length(fnames)
            par.(fnames{in})=parSource.(fnames{in});
        end
    catch
    end
end