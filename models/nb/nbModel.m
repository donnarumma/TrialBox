function   [data_out,out]=nbModel(data_in,par)
% function [data_out,out]=nbModel(dataTrials,par)

execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ', mfilename); end

labs        = [data_in.trialType]';
in_name     = par.InField;
data_3d     = cat(3,data_in.(in_name));
data_3d     = permute(data_3d,[3,1,2]);
data_2d     = reshape(data_3d,size(data_3d,1),size(data_3d,2)*size(data_3d,3));

cvp         = cvpartition(labs,'kfold',par.kfold,'Stratify',true);

mdl = fitcnb(data_2d, labs,...
    'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions', ...
    struct('Optimizer','randomsearch','CVPartition',cvp,'MaxObjectiveEvaluations',par.numIterations, 'AcquisitionFunctionName','expected-improvement-plus','ShowPlots',false,'Verbose',0));

data_out = data_in;
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
%% output save
out.mdl                     = mdl;
