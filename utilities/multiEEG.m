function data_out = multiEEG(data_trials,par)


Infield = par.Infield;

counter = 1;

data_out = struct();

for i = 1:length(data_trials)

    field_names = fieldnames(data_trials);

    current_eeg = data_trials(i).(Infield);

    num_slices = size(current_eeg, 1);

    for j = 1:num_slices
        for k = 1:length(field_names)
            field_name = field_names{k};
            if strcmp(field_name, Infield)
                data_out(counter,1).(field_name) = squeeze(current_eeg(j, :, :));
            else
                data_out(counter,1).(field_name) = data_trials(i).(field_name);
            end
        end
        data_out(counter,1).nEpoch = j;
        data_out(counter,1).OlDtrialId = data_out(counter).trialId;
        data_out(counter,1).trialId = counter;
        counter = counter + 1;
    end
end