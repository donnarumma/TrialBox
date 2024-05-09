function   [data_dir]=getDataSound()
% function [data_dir]=getDataSound()

[~, name] = system('cat /etc/os-release | grep -w NAME');
    
if contains(name, 'Ubuntu')
    data_dir = '~/DATA/';
else
    data_dir = 'D:\Articolo_Sound_Enrico\New_Extraction_files\';
end
end