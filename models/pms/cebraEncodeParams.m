function   params=cebraEncodeParams()
% function params=cebraEncodeParams()

params.exec                 = true;
params.xfld                 = 'time';
params.model_filename       = 'cebra_model.pkl';
params.group_field          = 'data';
params.manifold_field       = 'manifold';
params.manifold_filename    = 'manifold.hd5';
params.neural_field         = 'neural';
params.neural_filename      = 'neural.hd5';
params.script_filename      = 'cebraEncode.py';  % script to be executed in python (full path)
params.script_rundir        = './';             % directory where script expects inputs
params.encodeParams_filename= 'cebraEncodeParams.hd5';

