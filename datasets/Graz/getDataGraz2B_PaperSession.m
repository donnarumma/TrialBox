function   [data_dir]=getDataGraz2B_PaperSession()
% function [data_dir]=getDataGraz2B_PaperSession()

[~, name] = system('cat /etc/os-release | grep -w NAME');
    
if contains(name, 'Ubuntu')
    data_dir = '~/DATA/GRAZ/BCICIV_2b_mat/';
else
    data_dir = 'D:\D_Ausilio EEG\Graz_2b_paperSession';
end
end