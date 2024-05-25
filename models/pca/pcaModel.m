function   [data_trials, out]=pcaModel(data_trials,par)
% function [data_trials, out]=pcaModel(data_trials,par)

execinfo    =par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

InField     =par.InField;

Xdata      =[data_trials.(InField)]; % original data are nChannels x nTimes*nTrials 

% [dataTrials.(fn)] same as cat(2,dataTrials.(fn)) -> nChannels x nTimes*nTrials
% different from cat(3,dataTrials.(fn)) that put trials on the third dimension,
% i.e. nChannels x nTimes x nTrials 

% nChannels x nTimes*nTrials 
% -> this arrangement in pcaModel ignore that thata are collected in
% different time sequences and consideres each time step an independent snapshot 

% 
% nTimes*nTrials x nChannels, (set to perform matlab pca function - repetitions x variables)
% [Wpca, Zpca, ~, ~, explained, mu]
% do pca on the transpose
% Wpca:      nChannels x nComponents
% Zpca: nTimes*nTrials x nComponents
[Wpca, ~, ~, ~, explained] = pca(Xdata'); 
% Xdata Number of repetitions x number of variables (e.g. neurons). Thus mu=mean(Xdata,1);
% Zpca vectors the reduction space. Wpca coefficents
% data can be reconstructed as:
%  Zpca_rec        =     Xdata * Wpca
%  X_data_rec      = Zpca_data * Wpca' + mu
% inv(Wpca(mi,:)' * Wpca(mi,:)) * Wpca(mi,:)' * bsxfun(@minus, Y(mi,:), d(mi));
mu                            = mean(Xdata,2);
explain=cumsum(explained);
if par.numComponents<1
    whoexplain=find(explain>par.perc);
    par.numComponents=whoexplain(1);
end
numComponents=par.numComponents;
fprintf('Components: %g - Explained Variance on Data %g%%\n',numComponents,explain(numComponents));

% Zpca=Zpca(:,1:numComponents)';                          % numComponents x all times 
% Zpca=reshape(Zpca,numComponents,T,length(data_train));  % K x T x N trials in train
Wpca         = Wpca(:, 1:numComponents);                          % N x K
% if ~isempty(data_test)
%     Zpca_test   = (X_test-repmat(mu,1,size(X_test,2)))' * Wpca;                                         % T_test x K
%     Zpca_test   = reshape(Zpca_test',numComponents,T,length(data_test));  % K x T x N trials in test
% end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
out.Wpca          = Wpca;
out.mu            = mu;
out.explained     = explained;
out.numComponents = numComponents;