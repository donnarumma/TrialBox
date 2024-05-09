function   [data_dir]=getDataFITTSepotarget()
% function [data_dir]=getDataFITTSepotarget()

[~, name] = system('cat /etc/os-release | grep -w NAME');
    
if contains(name, 'Ubuntu')
    data_dir = '';
else
    data_dir = 'D:\D_Ausilio EEG\epotarget';
end
end