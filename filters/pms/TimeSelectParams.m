function   par = TimeSelectParams(parSource)
% function par = TimeSelectParams(parSource)
par.exec    = 1;
par.t1      = 0;
par.t2      = 0;
par.deltaT  = 0;
par.dt      = 0;

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