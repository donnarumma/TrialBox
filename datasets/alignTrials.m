function   [data, out] = alignTrials(data,params)
% function [data, out] = alignTrials(data,params)

execinfo=params.exec;
if ~isempty(execinfo); t=tic; end

id_event      = params.events; %TargetON
InField     = params.InField;
TimeField   = (['time' InField]);

events_time  =  data.events;
events_names = fieldnames(events_time);

for iTr = 1:length(data)
    offset = data(iTr).events.(id_event);
    data(iTr).(TimeField) = data(iTr).(TimeField) - data(iTr).events.(id_event);
    for n_ev = 1:length(events_names)
        data(iTr).events.(events_names{n_ev}) = data(iTr).events.(events_names{n_ev}) - offset;
    end
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end
