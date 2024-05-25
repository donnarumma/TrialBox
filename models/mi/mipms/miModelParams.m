function   par = miModelParams(parSource)
% function par = miModelParams(parSource)
par.exec            = true;
par.m               = 2;
par.k               = 4;
par.InField         = 'y';
par.OutField        = 'y';
par.csp_all         = false;
par.xfld        = 'time';
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end