function   [data_dir]=getDataGraz2B()
% function [data_dir]=getDataGraz()

[~, name] = system('cat /etc/os-release | grep -w NAME');
    
if contains(name, 'Ubuntu')
    data_dir = '~/DATA/GRAZ/BCICIV_2b_mat/';
else
    data_dir = 'D:\D_Ausilio EEG\Graz_2b_Extracted_data';
end
end