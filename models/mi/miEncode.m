function [data,out] = miEncode(data,par)

execinfo=par.exec;
if ~isempty(execinfo); t=tic; end

InField             = par.InField;
OutField            = par.OutField;
IndMI               = par.IndMI;
xfld                = par.xfld;


data_project        = cat(1,data.(InField));

NTrials             = size(data_project,1);

class = unique([data.trialType]');
n_class = length(class);
data_class = cell(n_class,1);
for ncl=1:n_class
    start_index = (ncl - 1) * size(data_project,2)/n_class + 1;
    end_index = ncl * size(data_project,2)/n_class;
    data_class{ncl} = data_project(:,start_index:end_index);
end

for ncl = 1:n_class
    data_class{ncl} = data_class{ncl}(:,IndMI{ncl});
end

data_project = horzcat(data_class{:});

for iTr=1:NTrials
    data(iTr).(OutField) = data_project(iTr,:);
    data(iTr).([xfld OutField])=data(iTr).([xfld InField])(end);
end

out.V = data_class;

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end