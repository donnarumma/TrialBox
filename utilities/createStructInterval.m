function [ResultKappa,ResultAcc] = createStructInterval(res,par)

subj = par.subj;
file = par.file;
% train_interval = [tr1 tr2];
train_name = par.train_name;
tr1 = par.train_tr1;
tr2 = par.train_tr2;
% test_interval = [ts1 ts2];
test_name = par.test_name;
ts1 = par.test_ts1;
ts2 = par.test_ts2;
class = par.class;
method = par.method;

try m = par.m;
catch
    m=0;
end

% KappaValue
ResultKappa.date = datetime('now');
ResultKappa.method = method;
ResultKappa.subj = subj;
ResultKappa.file = file;
ResultKappa.train_name = train_name;
ResultKappa.train_start = tr1;
ResultKappa.train_stop = tr2;
ResultKappa.train_Kappa_Prod = res.train.kappaValueProd;
ResultKappa.train_Kappa_Sum = res.train.kappaValueSum;
ResultKappa.train_Kappa_Max = res.train.kappaValueMax;
ResultKappa.test_name = test_name;
ResultKappa.test_start = ts1;
ResultKappa.test_stop = ts2;
ResultKappa.test_Kappa_Prod = res.test.kappaValueProd;
ResultKappa.test_Kappa_Sum = res.test.kappaValueSum;
ResultKappa.test_Kappa_Max = res.test.kappaValueMax;
ResultKappa.m = m;
for ncl=1:length(class)
    class_name = sprintf('Class%d',ncl);
    ResultKappa.(class_name) = class{ncl};
end
ResultKappa.irng = par.irng;

% Accuracy
ResultAcc.date = datetime('now');
ResultAcc_class.Method = method;
ResultAcc.subj = subj;
ResultAcc.file = file;
ResultAcc.train_name = train_name;
ResultAcc.train_start = tr1;
ResultAcc.train_stop = tr2;
ResultAcc.train_Acc_Prod = res.train.AccuracyProd;
ResultAcc.train_Acc_Sum = res.train.AccuracySum;
ResultAcc.train_Acc_Max = res.train.AccuracyMax;
ResultAcc.test_name = test_name;
ResultAcc.test_start = ts1;
ResultAcc.test_stop = ts2;
ResultAcc.test_Acc_Prod = res.test.AccuracyProd;
ResultAcc.test_Acc_Sum = res.test.AccuracySum;
ResultAcc.test_Acc_Max = res.test.AccuracyMax;
ResultAcc.m = m;
for ncl=1:length(class)
    class_name = sprintf('Class%d',ncl);
    ResultAcc.(class_name) = class{ncl};
end
ResultAcc.irng = par.irng;
