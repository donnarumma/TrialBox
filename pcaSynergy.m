function   [dataTrials, out]=pcaSynergy(dataTrials,par)
% function [dataTrials, out]=pcaCompute(dataTrials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; end

InField       = par.InField;
OutField      = par.OutField;

xfld = 'time';
if length(size(dataTrials(1).(InField))) > 2
    dataTrain_3d = cat(4,dataTrials([dataTrials.train]).(InField));
    dataTrain_3d = permute(dataTrain_3d,[4,1,2,3]);
    dataTrain_3d  = reshape(dataTrain_3d,size(dataTrain_3d,1),size(dataTrain_3d,2)*size(dataTrain_3d,3)*size(dataTrain_3d,4)); % trials x variables*time
else
    dataTrain_3d  = cat(3,dataTrials([dataTrials.train]).(InField));                                      % variables x time x trials.
    dataTrain_3d  = permute(dataTrain_3d,[3,1,2]);                                                        % trials x variables x time
    dataTrain_3d  = reshape(dataTrain_3d,size(dataTrain_3d,1),size(dataTrain_3d,2)*size(dataTrain_3d,3)); % trials x variables*time
end
if length(size(dataTrials(1).(InField))) > 2
    dataTest_3d = cat(4,dataTrials([dataTrials.test]).(InField));
    dataTest_3d = permute(dataTest_3d,[4,1,2,3]);
    dataTest_3d  = reshape(dataTest_3d,size(dataTest_3d,1),size(dataTest_3d,2)*size(dataTest_3d,3)*size(dataTest_3d,4)); % trials x variables*time
else
    dataTest_3d   = cat(3,dataTrials([dataTrials.train]).(InField));                                      % variables x time x trials.
    dataTest_3d   = permute(dataTest_3d,[3,1,2]);                                                         % trials x variables x time
    dataTest_3d   = reshape(dataTest_3d,size(dataTest_3d,1),size(dataTest_3d,2)*size(dataTest_3d,3));     % trials x variables*time
end

% [Wpca, Zpca, ~, ~, explained, mu]
[Wpca_train, Zpca_train, ~, ~, explained] = pca(dataTrain_3d); % data_3d Number of trials x variables*time. Thus mu=mean(Xdata,1);
% Zpca vectors the reduction space. Wpca coefficents
% data can be reconstructed as:
%  Zpca_rec        =     Xdata * Wpca
%  X_data_rec      = Zpca_data * Wpca' + mu
% inv(Wpca(mi,:)' * Wpca(mi,:)) * Wpca(mi,:)' * bsxfun(@minus, Y(mi,:), d(mi));

explain  = cumsum(explained);
if par.numComponents<1
    whoexplain=find(explain>par.perc);
    par.numComponents=whoexplain(1);
end
numComponents=par.numComponents;
fprintf('Components: %g - Explained Variance on Data %g%%\n',numComponents,explain(numComponents));

Wpca_train=Wpca_train(:, 1:numComponents);      % nvar          x numcomponents
Zpca_train=Zpca_train(:,1:numComponents);       % numComponents x all times
indsTrain=find([dataTrials.train]);
for iT = 1:length(indsTrain)
    % dataTrials(indsTrain(iT)).(par.OutField) = Zpca_train(iT,:)';
    dataTrials(indsTrain(iT)).(par.OutField) = Zpca_train(iT,:);
    dataTrials(iT).([xfld OutField])=dataTrials(iT).([xfld InField])(end);
end
indsTest=find([dataTrials.test]);
if ~isempty(indsTest)
    Zpca_test=dataTest_3d*Wpca_train;
    for iT = 1:length(indsTest)
        % dataTrials(indsTest(iT)).(par.OutField) = Zpca_test(iT,:)';
        dataTrials(indsTest(iT)).(par.OutField) = Zpca_test(iT,:);
    end
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end
out.W             =Wpca_train;
out.explained     =explained;
out.numComponents =numComponents;