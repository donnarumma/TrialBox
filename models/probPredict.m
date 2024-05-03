function Results = probPredict(data_trials,par)

execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ', mfilename); end

InField                 = par.InField;
labs_train              = par.label_train;
labs_test               = par.label_test;

% TRAIN
Train_values = struct();
for nTrl=1:length(labs_train)
    for nInt=1:length(data_trials)
        Train_values(nTrl).Prob(:,nInt) = data_trials(nInt).(InField).train.prob(nTrl,:)';
    end
end

pTrain_pred = struct();
labsTrain_pred = struct();
for ntr=1:length(labs_train)
    pTrain_prod          = prod(Train_values(ntr).Prob,2);
    pTrain_final_prod    = pTrain_prod/sum(pTrain_prod);
    pTrain_sum           = sum(Train_values(ntr).Prob,2);
    pTrain_final_sum     = pTrain_sum/sum(pTrain_sum);
    pTrain_max           = max(Train_values(ntr).Prob,[],2);
    pTrain_final_max     = pTrain_max/sum(pTrain_max);
    % Prediction probability and labels
    [pTrain_pred(ntr).prod, labsTrain_pred(ntr).prod]  = max(pTrain_final_prod);
    [pTrain_pred(ntr).sum, labsTrain_pred(ntr).sum]    = max(pTrain_final_sum);
    [pTrain_pred(ntr).max, labsTrain_pred(ntr).max]    = max(pTrain_final_max);
end
% Product Probability
Results.train.AccuracyProd                  = sum([labsTrain_pred.prod]' == labs_train)/length(labs_train)*100;
CmatrxixTrain_Prod                          = confusionmat(labs_train, [labsTrain_pred.prod]');
Results.train.kappaValueProd                = kappaModel(CmatrxixTrain_Prod);
Results.train.AccuracyProd_class(1,:)       = accuracy4classes(labs_train,[labsTrain_pred.prod]');

% Sum Probability
Results.train.AccuracySum                   = sum([labsTrain_pred.sum]' == labs_train)/length(labs_train)*100;
CmatrxixTrain_Sum                           = confusionmat(labs_train,[labsTrain_pred.sum]');
Results.train.kappaValueSum                 = kappaModel(CmatrxixTrain_Sum);
Results.train.AccuracySum_class(1,:)        = accuracy4classes(labs_train,[labsTrain_pred.sum]');

% Max Probability
Results.train.AccuracyMax                   = sum([labsTrain_pred.max]' == labs_train)/length(labs_train)*100;
CmatrxixTrain_Max                           = confusionmat(labs_train, [labsTrain_pred.max]');
Results.train.kappaValueMax                 = kappaModel(CmatrxixTrain_Max);
Results.train.AccuracyMax_class(1,:)        = accuracy4classes(labs_train,[labsTrain_pred.max]');


% TEST
Test_values = struct();
for nTrl=1:length(labs_test)
    for nInt=1:length(data_trials)
        Test_values(nTrl).Prob(:,nInt) = data_trials(nInt).(InField).test.prob(nTrl,:)';
    end
end

pTest_pred = struct();
labsTest_pred = struct();
for ntr=1:length(labs_test)
    pTest_prod          = prod(Test_values(ntr).Prob,2);
    pTest_final_prod    = pTest_prod/sum(pTest_prod);
    pTest_sum           = sum(Test_values(ntr).Prob,2);
    pTest_final_sum     = pTest_sum/sum(pTest_sum);
    pTest_max           = max(Test_values(ntr).Prob,[],2);
    pTest_final_max     = pTest_max/sum(pTest_max);
    % Prediction probability and labels
    [pTest_pred(ntr).prod, labsTest_pred(ntr).prod]  = max(pTest_final_prod);
    [pTest_pred(ntr).sum, labsTest_pred(ntr).sum]    = max(pTest_final_sum);
    [pTest_pred(ntr).max, labsTest_pred(ntr).max]    = max(pTest_final_max);
end

% Product Probability
Results.test.AccuracyProd                   = sum([labsTest_pred.prod]' == labs_test)/length(labs_test)*100;
CmatrxixTest_Prod                           = confusionmat(labs_test,[labsTest_pred.prod]');
Results.test.kappaValueProd                 = kappaModel(CmatrxixTest_Prod);
Results.test.AccuracyProd_class(1,:)        = accuracy4classes(labs_test,[labsTest_pred.prod]');

% Sum Probability
Results.test.AccuracySum                    = sum([labsTest_pred.sum]' == labs_test)/length(labs_test)*100;
CmatrxixTest_Sum                            = confusionmat(labs_test,[labsTest_pred.sum]' );
Results.test.kappaValueSum                  = kappaModel(CmatrxixTest_Sum);
Results.test.AccuracySum_class(1,:)         = accuracy4classes(labs_test,[labsTest_pred.sum]');

% Max Probability
Results.test.AccuracyMax                    = sum([labsTest_pred.max]' == labs_test)/length(labs_test)*100;
CmatrxixTest_Max                            = confusionmat(labs_test, [labsTest_pred.max]');
Results.test.kappaValueMax                  = kappaModel(CmatrxixTest_Max);
Results.test.AccuracyMax_class(1,:)         = accuracy4classes(labs_test,[labsTest_pred.max]');

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
