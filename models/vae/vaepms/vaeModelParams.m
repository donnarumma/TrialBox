function   par = vaeModelParams(parSource)
% function par = vaeModelParams(parSource)
par.exec                = true;
par.InField             = 'y';
par.numEpochs           = 150;
par.miniBatchSize       = 128;
par.learnRate           = 1e-3;
par.numLatentChannels   = 32;
par.kernsize            = 3;
par.projsize            = 7;                
par.MBnumOutputs        = 1;
par.validPercent        = 0;
par.validationFrequency = 1;
par.graphplot           = true;
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end