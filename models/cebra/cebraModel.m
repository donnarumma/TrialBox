function   [data_trials,out] = cebraModel(data_trials,params)
% function [data_trials,out] = cebraModel(data_trials,params)
execinfo            = params.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

% nTrials             = length(data_trials);
script_filename     = params.script_filename;
script_rundir       = params.script_rundir;
group_field         = params.group_field;
neural_field        = params.neural_field;
behavior_field      = params.behavior_field;
neural              = [data_trials.(neural_field)];    % nChannels x nTimes*nTrials
behavior            = [data_trials.(behavior_field)];  % nFeatures x nTimes*nTrials

% https://stackoverflow.com/questions/21624653/python-created-hdf5-dataset-transposed-in-matlab
% note: same hd5 file results transposed 
% in matlab: Channels x nTimes
% in python:   nTimes x nChannels 
neural_filename     = [script_rundir params.neural_filename];
fprintf('Saving %s\n',neural_filename)
params.neural_filename= neural_filename;
h5create(neural_filename,['/' group_field '/' neural_field], size(neural));
h5write (neural_filename,['/' group_field '/' neural_field], neural);

% in matlab: nFeatures x nTimes
% in python:    nTimes x nFeatures 
behavior_filename       = [script_rundir params.behavior_filename];
fprintf('Saving %s\n',behavior_filename)
params.behavior_filename= behavior_filename;
h5create(behavior_filename,['/' group_field '/', behavior_field], size(behavior));
h5write (behavior_filename,['/' group_field '/', behavior_field], behavior);

modelParams_filename                = [params.script_rundir params.modelParams_filename];
fnames                              = fieldnames(params);
if ~isfile(modelParams_filename)
    h5create(modelParams_filename, '/dummy', [1 1]);  % Create a dummy dataset
end
fid = H5F.open(modelParams_filename,'H5F_ACC_RDWR','H5P_DEFAULT');
H5L.delete(fid,'/dummy','H5P_DEFAULT');
for ifld = 1:length(fnames)
    fattr   = params.(fnames{ifld});
    if islogical(fattr)
        if fattr
            fattr='true';
        else
            fattr='false';
        end
    end
    h5writeatt(modelParams_filename,'/',fnames{ifld},fattr);
end
H5F.close(fid);
% h5disp(modelParams_filename)
%% run fit in python
modelCommand        = ['python ' script_filename ' ' modelParams_filename];
fprintf('Executing %s\n',modelCommand)
[~,message]         = system(modelCommand);
disp(message);

out.model_filename  = params.model_filename;

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
end
