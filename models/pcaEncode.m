function   [data_trials, out]=pcaEncode(data_trials,par)
% function [data_trials, out]=pcaEncode(data_trials,par)
% pca executed per class
% inspired by demixed pca https://elifesciences.org/articles/10989 
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
nTrials     = length(data_trials);

InField     = par.InField;
OutField    = par.OutField;

[~,nTimes]       = size(data_trials(1).(InField));
X_data      = [data_trials.(InField)];          % data are nChannels x nTimes * nTrials 
X_data      = X_data - repmat(par.mu,1,size(X_data,2));

% [dataTrials.(InField)] same as cat(2,dataTrials.(fn)) 
% different from cat(3,dataTrials.(fn)) that put trials on the third dimension
% transpose to set X_data nTimes x nChannels, in order to perform a pca on each time step

% do pca on the transpose (nTimes*nTrials x nChannels)
% Wpca: nChannels x nComponents
% Zpca: nTimes x nComponents
% [Wpca, ~, ~, ~, explained] = pca(X_data');
% Xdata Number of repetitions x number of variables (e.g. neurons). Thus mu=mean(Xdata,1);
% Zpca vectors the reduction space. Wpca coefficents
% data can be reconstructed as:
%  Zpca_rec        =     X_data' * Wpca  % TxN x NxK 
%  X_data_rec      = Zpca_data  * Wpca' + mu
% inv(Wpca(mi,:)' * Wpca(mi,:)) * Wpca(mi,:)' * bsxfun(@minus, Y(mi,:), d(mi));
explained              =par.explained;
Wpca                   =par.Wpca;
[~,nComponents]        =size(Wpca);
explain                =cumsum(explained);


fprintf('Components: %g - Explained Variance on Data %g%%\n',nComponents,explain(nComponents));

Wpca=Wpca(:, 1:nComponents);                           % nChannels x nComponents

Zpca    = X_data' * Wpca;                              % nTimes*nTrials x nComponents
Zpca    = reshape(Zpca',nComponents,nTimes,nTrials);   % nComponents x nTimes x nTrials

xfld    = 'time';
for iTrial=1:nTrials
    data_trials(iTrial).(OutField)          = Zpca(:,:,iTrial);
    data_trials(iTrial).([xfld OutField])   = data_trials(iTrial).([xfld InField]);
end


if ~isempty(execinfo); out.exectime=toc(t); fprintf('Time Elapsed: %.2f s\n',out.exectime); end
out.explained     =explained;
out.numComponents =nComponents;