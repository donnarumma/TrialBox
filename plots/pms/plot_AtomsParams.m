function   par = plot_AtomsParams(parSource)
% function par = plot_AtomsParams(parSource)
par.exec        = true;
par.hfig        = [];
par.cmaps       = [];
par.InField     = 'dict';
par.nRows       = 10;
par.xfld        = 'time';
par.Aname       = 'A';
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end