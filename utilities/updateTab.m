
function updated_table = updateTab(data,params)

% function updated_table = updateTab(data,params)
% This function can update an existing excel file with new results

dir = params.dir;
name = params.name;
sheetnames = params.sheetnames;
% % CSV
% filename = [dir,'\',name,'.csv'];
% XLSX
filename = [dir,'\',name,'.xlsx'];

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

updated_table = [existing_table; new_table];

% Save Table with new Results
% XLSX
% writetable(updated_table, filename, 'Sheet', sheetnames, 'WriteMode', 'append');
writetable(updated_table, filename, 'Sheet', sheetnames, 'WriteMode', 'overwritesheet');

% % CSV
% writetable(updated_table, filename, 'WriteMode', 'append');