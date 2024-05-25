function   par = faModelParams(parSource)
% function par = faModelParams(parSource)
par.exec            = true;
par.InField         = 'y';
par.OutField        = 'y';
par.tol             = 1e-8; 
par.nIters          = 1e8;   % number of EM steps
par.minVarFrac      = 0.01;
par.numComponents   = 3;
% choose between 'fa'  -> factor analysis
% choose between 'ppa' -> probabilistic pca
par.type            = 'fa';  
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end