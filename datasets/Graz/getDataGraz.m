function   [data_dir]=getDataGraz()
% function [data_dir]=getDataGraz()

[~, name] = system('cat /etc/os-release | grep -w NAME');
    
if contains(name, 'Ubuntu')
    data_dir = '~/DATA/GRAZ/BCICIV_2a_mat/';
else
    data_dir = 'D:\D_Ausilio EEG\EEG_FITTS\Data_Graz_Extracted\Extracted_Data\Data_Paper_Interval';
end
end