function   [data_dir]=getDataAudiobook()
% function [data_dir]=getDataAudiobook()

[~, name] = system('cat /etc/os-release | grep -w NAME');
    
if contains(name, 'Ubuntu')
    data_dir = '/home/frosolone/HILARIOUS_Sound/Audiobook';
else
    data_dir = 'E:\Hilarious_Dubbing\Audiobook\';
end
end