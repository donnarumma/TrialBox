function   par = vaeDecodeParams(parSource)
% function par = vaeDecodeParams(parSource)
par.exec                = true;
par.InField             = 'y';
par.OutField            = 'y';
par.xfld                = 'time';
par.miniBatchSize       = 128;
par.MBnumOutputs        = 1;
par.netD                = [];
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end