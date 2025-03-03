
function updated_table = updateTab(data,params)

% function updated_table = updateTab(data,params)
% This function can update an existing excel file with new results

dir = params.dir;
name = params.name;
sheetnames = params.sheetnames;

[~, namedir] = system('cat /etc/os-release | grep -w NAME');
if contains(namedir, 'Ubuntu')
    sep = '/';
else
    sep ='\';
end

% % CSV
% filename = [dir,'\',name,'.csv'];
% XLSX
filename = [dir,sep,name,'.xlsx'];

try
    existing_table = readtable(filename,'Sheet', sheetnames);
catch
    T = table();
    % % CSV
    % writetable(T, filename)
    % existing_table = readtable(filename);
    % XLSX
    writetable(T, filename, 'Sheet',sheetnames);
    existing_table = readtable(filename,'Sheet', sheetnames);
end

new_table = struct2table(data);
existing_table = table2struct(existing_table);

if isfield(existing_table, 'train_start') || isfield(existing_table, 'ProbabilityMean')

    if ~ischar(existing_table(1).train_start)
        for i=1:length(existing_table)
            existing_table(i).train_start =  sprintf('%s  ', num2str(existing_table(i).train_start));
            existing_table(i).train_stop =  sprintf('%s  ', num2str(existing_table(i).train_stop));
            existing_table(i).test_start =  sprintf('%s  ', num2str(existing_table(i).test_start));
            existing_table(i).test_stop =  sprintf('%s  ', num2str(existing_table(i).test_stop));
            existing_table(i).ProbabilityMean  =  sprintf('%s  ', num2str(existing_table(i).ProbabilityMean));

        end
    end
end
existing_table = struct2table(existing_table);


% % Trova le colonne mancanti nella nuova tabella rispetto alla tabella esistente
% A = existing_table.Properties.VariableNames;
% B =  new_table.Properties.VariableNames;
% missing_columns = setdiff(A,B);
%
% % Aggiungi le colonne mancanti alla nuova tabella e impostale a NaN
% for i = 1:numel(missing_columns)
%     new_table.(missing_columns{i}) = NaN(size(new_table, 1), 1);
% end
try
    updated_table = [existing_table; new_table];
catch
    new_table = table2struct(new_table);
    existing_table = table2struct(existing_table);
            for i=1:length(new_table)
                new_table(i).train_start =  sprintf('%s', num2str(existing_table(i).train_start));
                new_table(i).train_stop =  sprintf('%s', num2str(existing_table(i).train_stop));
                new_table(i).test_start =  sprintf('%s', num2str(existing_table(i).test_start));
                new_table(i).test_stop =  sprintf('%s', num2str(existing_table(i).test_stop));
            end
    new_table = struct2table(new_table);
    existing_table = struct2table(existing_table);

    updated_table = [existing_table; new_table];
end
% Save Table with new Results
% XLSX
% writetable(updated_table, filename, 'Sheet', sheetnames, 'WriteMode', 'append');
% writetable(updated_table, filename, 'Sheet', sheetnames, 'WriteMode', 'overwritesheet');
writetable(updated_table, filename, 'Sheet', sheetnames, 'WriteMode', 'overwritesheet', 'WriteVariableNames', true, 'WriteRowNames', true, 'FileType', 'spreadsheet');
% % CSV
% writetable(updated_table, filename, 'WriteMode', 'append');