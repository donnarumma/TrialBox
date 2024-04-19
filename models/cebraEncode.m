function   [data_trials,out] = cebraEncode(data_trials, params)
% function [data_trials,out] = cebraEncode(data_trials, par)
execinfo            = params.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

script_filename         = params.script_filename;
script_rundir           = params.script_rundir;
group_field             = params.group_field;
neural_field            = params.neural_field;
manifold_field          = params.manifold_field;
xfld                    = params.xfld;

neural                  = [data_trials.(neural_field)];    % nChannels x nTimes*nTrials
% create file of inputs
neural_filename         = [script_rundir params.neural_filename];
fprintf('Saving %s\n',neural_filename)
params.neural_filename  = neural_filename;
h5create(neural_filename,['/' group_field '/' neural_field], size(neural));
h5write (neural_filename,['/' group_field '/' neural_field], neural);

% change put path to manifold file in script_rundir
manifold_filename       = [script_rundir params.manifold_filename];
params.manifold_filename= manifold_filename;

% create files of params
encodeParams_filename               = [params.script_rundir params.encodeParams_filename];
fnames                              = fieldnames(params);
if ~isfile(encodeParams_filename)
    h5create(encodeParams_filename, '/dummy', [1 1]);  % Create a dummy dataset
end
fid = H5F.open(encodeParams_filename,'H5F_ACC_RDWR','H5P_DEFAULT');
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
    h5writeatt(encodeParams_filename,'/',fnames{ifld},fattr);
end
H5F.close(fid);

% get trial size for each trial
nTrials                                     = length(data_trials);
trial_size                                  = nan(nTrials,1);
for iTrial=1:nTrials
    trial_size(iTrial)                      = size(data_trials(iTrial).(neural_field),2);
end

% run transform in python
encodeCommand       = ['python ' script_filename ' ' encodeParams_filename];
fprintf('Executing %s\n',encodeCommand)
[~,message]         = system(encodeCommand);
disp(message);

% load output result
fprintf('Loading %s\n',manifold_filename)
manifold                                    = h5read(manifold_filename, ['/' group_field '/' manifold_field]);

% put manifold in data_trials
bStart                                      = 1;
bEnd                                        = 0;
for iTrial = 1:nTrials
    bSize                                   = trial_size(iTrial);   % size of the current trial 
    bEnd                                    = bEnd + bSize;         % final index of the trial
    data_trials(iTrial).(manifold_field)    = manifold(:,bStart:bEnd); 
    bStart                                  = bEnd + 1;             % initial index of next trial

    data_trials(iTrial).([xfld neural_field]) = data_trials.([xfld neural_field]);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
end

