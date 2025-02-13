function   [data_dir]=getDataAudiobook()
% function [data_dir]=getDataAudiobook()

[~, name] = system('cat /etc/os-release | grep -w NAME');
    
if contains(name, 'Ubuntu')
    data_dir = '/home/frosolone/HILARIOUS_Sound/Audiobook';
else
    data_dir = 'D:\main_scriptNSA\Older_and_Proof\Sound_pipeline\11_06\ERP_analysis\Audiobook';
end
end