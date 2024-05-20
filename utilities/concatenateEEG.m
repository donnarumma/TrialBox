function data_out = concatenateEEG(data_trials1,data_trials2,par)

InField     = par.InField;
fsample     = par.fsample;
xfld        = 'time';
data_out = data_trials1;
for iTr=1:length(data_trials1)
    data_out(iTr).(InField) = cat(2,data_trials1(iTr).(InField),data_trials2(iTr).(InField));
    data_out(iTr,1).(['time' InField])   = linspace(0,size(data_out(iTr).(InField),2)/fsample,size(data_out(iTr).(InField),2));
end
