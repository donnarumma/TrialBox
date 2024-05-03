function     [data_trials, out] = remapTypes(data_trials,par)
% function   [data_trials, out] = remapTypes(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

trialTypes          = [data_trials.trialType];
trialNames          = {data_trials.trialName};
[types,idx_types]   = unique(trialTypes);
names               = trialNames(idx_types);         
if isempty(par.selection)
    sel_types       = num2cell(types);
else
    sel_types       = par.selection;
end
if ~isempty(par.names) & length(sel_types)==length(par.names)
    sel_names         = par.names;
else
    sel_names         = names([sel_types{:}]);
end

retain                = false(1,length(trialTypes));
for new_type = 1:length(sel_types)
    new_name = sel_names{new_type};
    lgc_sel = ismember(trialTypes,sel_types{new_type});
    retain = retain | lgc_sel;
    indx = find(lgc_sel);
    for it=indx
        data_trials(it).trialType  =new_type;
        data_trials(it).trialName  =new_name;
    end
end
data_trials=data_trials(retain);
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end