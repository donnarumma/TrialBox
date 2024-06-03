function   par = reconstructionErrorsParams(parSource)
% function par = reconstructionErrorsParams(parSource)
par.exec                = true;
par.FieldReconstructed  = 'y_rec';
par.FieldTarget         = 'y';
    try
        fnames=fieldnames(parSource);
        for in=1:length(fnames)
            par.(fnames{in})=parSource.(fnames{in});
        end
    catch
    end
end