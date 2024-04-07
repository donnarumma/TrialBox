function   [data_trials, out]=pcaEncode(data_trials,par)
% function [data_trials, out]=pcaEncode(data_trials,par)
% pca executed per class
% inspired by demixed pca https://elifesciences.org/articles/10989 
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
nTrials     = length(data_trials);

InField     = par.InField;
OutField    = par.OutField;

% data_train  =data_trials([data_trials.train]);
% data_test   =data_trials([data_trials.test]); % original data are N variables x Time T. 

[~,T]       = size(data_trials(1).(InField));
X_data      = [data_trials.(InField)];          % data are nChannels x nTimes * nTrials 
X_data      = X_data - repmat(par.mu,1,size(X_data,2));

% [dataTrials.(fn)] same as cat(2,dataTrials.(fn)) different from cat(3,dataTrials.(fn)) that put trials on the third dimension
% transpose to set Xdata T x variablex, in order to perform a pca on each time step
% 
% labels_train=[data_train.trialType];
% classes     =unique(labels_train);
% Nclasses    =length(classes);
% X_train_c    =nan(Nvar,T,Nclasses);

% X_train_c            =reshape(X_train_c,Nvar,T*Nclasses); % Nvar x T*Nclasses
% % [Wpca, Zpca, ~, ~, explained, mu]
% [Wpca, ~, ~, ~, explained] = pca(X_train_c'); % do pca on the transpose (TxN). Wpca: N x K. Zpca: TxK
% Xdata Number of repetitions x number of variables (e.g. neurons). Thus mu=mean(Xdata,1);
% Zpca vectors the reduction space. Wpca coefficents
% data can be reconstructed as:
%  Zpca_rec        =     Xdata' * Wpca    % TxN x NxK 
%  X_data_rec      = Zpca_data  * Wpca' + mu
% inv(Wpca(mi,:)' * Wpca(mi,:)) * Wpca(mi,:)' * bsxfun(@minus, Y(mi,:), d(mi));
explained              =par.explained;
Wpca                   =par.Wpca;
[~,nComponents]        =size(Wpca);
explain                =cumsum(explained);


fprintf('Components: %g - Explained Variance on Data %g%%\n',nComponents,explain(nComponents));

Wpca=Wpca(:, 1:nComponents);                                            % N x K

Zpca    = X_data' * Wpca;                                           % T_train x K
Zpca    = reshape(Zpca',nComponents,T,nTrials);   % K x T x N trials in test

xfld    = 'time';
for iTrial=1:nTrials
    data_trials(iTrial).(OutField)          = Zpca(:,:,iTrial);
    data_trials(iTrial).([xfld OutField])   = data_trials(iTrial).([xfld InField]);
end


if ~isempty(execinfo); out.exectime=toc(t); fprintf('Time Elapsed: %.2f s\n',out.exectime); end
% out.W             =Wpca;
% out.mu            =mean(X_train_c,2);
out.explained     =explained;
out.numComponents =nComponents;