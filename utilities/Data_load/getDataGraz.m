function   [data_dir]=getDataGraz()
% function [data_dir]=getDataGraz()

[~, name] = system('cat /etc/os-release | grep -w NAME');
    
if contains(name, 'Ubuntu')
    data_dir = '~/DATA/GRAZ/BCICIV_2a_gdf/';
else
    data_dir = 'D:\D_Ausilio EEG\EEG_FITTS\Data_Graz_Extracted\Extracted_Data\PaperInterval';
end
end