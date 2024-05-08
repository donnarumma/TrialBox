function t = time_event(data_trials,idevent)

out = NaN(1,length(data_trials));
for iTr=1:length(data_trials)
    out(1,iTr) = data_trials(iTr).events.(idevent);
end

t = min(out);