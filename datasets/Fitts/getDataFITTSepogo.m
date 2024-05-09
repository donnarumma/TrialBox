function   [data_dir]=getDataFITTSepogo()
% function [data_dir]=getDataFITTSepogo()

[~, name] = system('cat /etc/os-release | grep -w NAME');
    
if contains(name, 'Ubuntu')
    data_dir = '';
else
    data_dir = 'D:\D_Ausilio EEG\epogo';
end
end