function   [data_trials,out] = cebraEncode(data_trials, par)
% function [data_trials,out] = cebraEncode(data_trials, par)
execinfo                                    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
model_file                                  = par.model_file;
InField                                     = par.InField;
OutField                                    = par.OutField;
xfld                                        = par.xfld; %'time';

% matrici vuote per mettere dati ratto
% data_trials=rat_data;
data_n                                      = [];
nTrials                                     = length(data_trials);
trial_size                                  = nan(nTrials,1);
for iTrial=1:nTrials
    data_n                                  = [data_n; data_trials(iTrial).(InField)'];
    % save trial dimensions
    trial_size(iTrial)                      = size(data_trials(iTrial).(InField),2);
end
script_file                                 = par.script_transform;
script_input_dir                            = par.script_input_dir;
script_output_dir                           = par.script_output_dir;
save([script_input_dir,'data_n.mat'], 'data_n');

% project data in python
command = sprintf('python "%s" "%s" "%s" "%s"', script_file, model_file, script_input_dir, script_output_dir);
fprintf('%s', command);  
[~, cmdout]                                 = system(command);

fprintf('%s', cmdout);  

out.model_transform                         = load(fullfile(script_output_dir, 'transf_data.mat'));
transf_data                                 = out.model_transform.transformed_data;

bStart                                      = 1;
bEnd                                        = 0;
for iTrial = 1:nTrials
    bSize                                   = trial_size(iTrial); % size of the current trial 
    bEnd                                    = bEnd + bSize; % Calcola l'indice finale per questo trial
    data_trials(iTrial).(OutField)          = transf_data(bStart:bEnd, :)'; 
    bStart                                  = bEnd + 1; %  indice iniziale per il prossimo trial

    data_trials(iTrial).([xfld OutField])   = data_trials.([xfld InField]);
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
end

