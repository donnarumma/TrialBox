function   [data, out] = tsneModel(data,par)
% function [data, out] = tsneModel(data,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

InField  =par.InField;
OutField =par.OutField;
xfld     =par.xfld;
data_all_3d         = cat(3,data.(InField));
X_data              = reshape(data_all_3d,size(data_all_3d,1)*size(data_all_3d,2),size(data_all_3d,3))';

nTrials             = length(data);
% t-sne
X_tsne              = tsne(X_data,'NumDimensions',par.NumDimensions);

for in = 1:nTrials
    data(in).(OutField) = X_tsne(in,:)';
    data(in).([xfld OutField])=data(in).([xfld InField])(end);
end

%% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end