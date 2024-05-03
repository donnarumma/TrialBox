function   par = remapTypesParams(parSource)
% function par = remapTypesParams(parSource)
par.exec      = true;
par.selection = [];
par.names     = [];
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