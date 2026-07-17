function params = getDefaultSPIKES_PLOTparams()
% function params = getDefaultSPIKES_PLOTparams()


%% --- DEFAULT SELECTION ---
params.monkey           = {'S','K'}; % selection Monkey "S", "K" or both

params.numComponents    = 8; % number of pca components
params.wd2D               = [3,4]; % pca component to 2D plot
params.wd3D               = [2,3,4]; % pca component to 3D plot


params.tplotStart       = -80; % plot start time
params.tplotEnd         = 200; % plot end time
params.keep             = [1:8;9:16;17:24]; % conditions 1:8 Act S - Obs K 
                                            % 9:16 Obs S - Act K
                                            % 17:24  Act S - Act K

params.binWidth         = 20;
params.kernSD           = 30;
%% SAVE PARAMS
params.dir_save         = 'SPIKES_PLOT';
