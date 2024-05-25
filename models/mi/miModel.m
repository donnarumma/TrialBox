function [data, out] = miModel(data,par)
% function Ind_final = miModel(data,par)
% Mutual Information (mi) executed per class
% inspired by MutualInformation algorithm https://ieeexplore.ieee.org/abstract/document/5332383
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s', mfilename); end

InField     = par.InField;

data_train = cat(1,data.(InField));
labs_train = [data.trialType]';
class = unique(labs_train);
n_class = length(class);
Ind_final = cell(n_class,1);
for ncl=1:n_class
    start_index = (ncl - 1) * size(data_train,2)/n_class + 1;
    end_index = ncl * size(data_train,2)/n_class;
    Ind_final{ncl} = MutualInformation(data_train(:,start_index:end_index,:),labs_train,ncl,par.m,par.k);
end

out.IndMI = Ind_final;
% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end