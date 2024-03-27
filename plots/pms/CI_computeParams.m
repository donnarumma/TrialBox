function   par = CI_computeParams(par)
% function par = CI_computeParams(par)

par.P    = 95;      % Confidence interval 68 (1 std)- 95 (2 std)- 99.7 (3 std)
par.N    = 10000;   % number of bootstrap samples
% params.opt=[1,0,0];  % BOOTSTRAP
% params.opt=[0,1,0];  % BOOTSTRAP, equivalent to previous option - slower and demonstrative
% params.opt=[0,0,1];  % mean and variance
par.opt  = [0,0,1]; % [1,0,0] perform bootstrap vs [0,0,1] without bootstrap
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end