function   [data_trials, out]=run_trials(data_trials,par)
% function [data_trials, execinfo]=run_trials(data_trials,par)
% execinfo = par.exec;
parfunname        = par.exec.funname;
parfun            = cellfun(@(X) (str2func(X)), parfunname,'uniformoutput',false);

for ifun=1:length(parfunname)
    parfield                = parfunname{ifun};
    par.(parfield).exec     = true;
    try
        income                  = par.(parfield).(income);
        outfun                  = par.(parfield).(outfun);
        outcome                 = par.(parfield).(outcome);
        par.(parfield).(income) = out.(outfun).(outcome);
    catch
    end
    % fprintf('Executing %s\n',parfield);
    % [data_trials, fout] = par.exec.fun{ifun}(data_trials,par.(parfield));
    [data_trials, fout]     = parfun{ifun}(data_trials,par.(parfield));
    % execinfo                = fout.execinfo;
    out.(parfield)          = fout;
     % execinfo.out
end