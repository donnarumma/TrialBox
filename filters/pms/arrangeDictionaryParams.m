function   par = arrangeDictionaryParams(parSource)
% function par = arrangeDictionaryParams(parSource)
par.exec        = true;
par.xfld        = 'time';
par.dictionary  = [];
par.nChannels   = 3;
par.nTimes      = 1;
par.OutField    = 'dict';
par.times       = [];
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end