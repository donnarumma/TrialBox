function   params=fccaModelParams()
% function params=fccaModelParams()

params.exec                 = true;
params.d                    = 3;
params.T                    = 4;
params.init                 = 'random_ortho';
params.n_init               = 50;
params.stride               = 11;
params.tol                  = 1e-6;
params.ortho_lambda         = 10.;
params.verbose              = 'false';
params.device               = 'cpu';
params.rng_or_seed          = 1;
% params.seed                 = 0;
% params.maxI                 = 9999;
params.group_field          = 'data';
params.neural_field         = 'neural';
params.model_field          = 'model';
params.neural_filename      = 'neural.hd5';
params.model_filename       = 'fcca_model.hd5';

params.script_filename      = 'fccaModel.py';   % script to be executed in python (full path)
params.script_rundir        = './';             % directory where script expects inputs
params.modelParams_filename = 'fccaModelParams.hd5';