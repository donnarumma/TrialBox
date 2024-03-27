function   [data_trials, out]=pcaModel(data_trials,par)
% function [data_trials, out]=pcaModel(dataTrials,par)

execinfo     =par.exec;
if ~isempty(execinfo); t=tic; end

InField     =par.InField;
% OutField    =par.OutField;


X_data      =[data_trials.(InField)]; % original data are N variables  x Time T. 
% X_test       =[data_test.(InField)];  % original data are N variables  x Time T. 

% [dataTrials.(fn)] same as cat(2,dataTrials.(fn)) different from cat(3,dataTrials.(fn)) that put trials on the third dimension
% transpose to set Xdata T x variablex, in order to perform a pca on each time step
% [Wpca, Zpca, ~, ~, explained, mu]
[Wpca, ~, ~, ~, explained] = pca(X_data'); % do pca on the transpose (TxN). Wpca: N x K. Zpca: TxK
% Xdata Number of repetitions x number of variables (e.g. neurons). Thus mu=mean(Xdata,1);
% Zpca vectors the reduction space. Wpca coefficents
% data can be reconstructed as:
%  Zpca_rec        =     Xdata * Wpca
%  X_data_rec      = Zpca_data * Wpca' + mu
% inv(Wpca(mi,:)' * Wpca(mi,:)) * Wpca(mi,:)' * bsxfun(@minus, Y(mi,:), d(mi));
mu                            = mean(X_data,2);
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

% allTrialsID =sort([data_trials.trialId]);
% data_trials =data_trials(allTrialsID);
% allTrialsID =[data_trials.trialId];
% % train
% trainIds    =[data_train.trialId];
% for iD=1:length(trainIds)
%     ipos=trainIds(iD)==allTrialsID;
%     data_trials(ipos).(OutField)=Zpca(:,:,iD);
%     data_trials(ipos).(['time' OutField])=data_trials(ipos).(['time' InField]);
% end
% % test 
% testIds    =[data_test.trialId];
% for iD=1:length(testIds)
%     ipos=testIds(iD)==allTrialsID;
%     data_trials(ipos).(OutField)=Zpca_test(:,:,iD);
%     data_trials(ipos).(['time' OutField])=data_trials(ipos).(['time' InField]);
% end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end
out.W             =Wpca;
out.mu            =mu;
out.explained     =explained;
out.numComponents =numComponents;