function   par = BootStrapDataParams(parSource)
% function par = BootStrapDataParams(parSource)
par.exec        = true;
par.ifclass     = true; % true  save output in data_trials     and original in out.data_trials
                        % false save output in out.data_trials and original data_trials
par.InField     = 'y';
par.OutField    = 'y';
par.P           = 95;      % Confidence interval 68 (1 std)- 95 (2 std)- 99.7 (3 std)
par.N           = 10000;   % number of bootstrap samples
% par.opt=[1,0,0];  % BOOTSTRAP
% par.opt=[0,1,0];  % BOOTSTRAP, equivalent to previous option - slower and demonstrative
% par.opt=[0,0,1];  % mean and variance
par.opt         = [0,0,1]; % [1,0,0] perform bootstrap vs [0,0,1] without bootstrap
par.SE          = 0; % 0 standard deviation, 1 standard error, bootstrap standard error 
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end