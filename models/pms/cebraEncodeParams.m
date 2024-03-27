function   params=cebraEncodeParams()
% function params=cebraEncodeParams()

params.exec             = true;
params.InField          = 'y';
params.OutField         = 'y';
params.script_transform = 'transf_py_rat.py';
params.script_output_dir= './';
params.script_input_dir = './';
params.xfld             = 'time';
