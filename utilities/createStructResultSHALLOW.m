function ResultAcc = createStructResultSHALLOW(Acc,par)

subj = par.subj;
file = par.file;
% train_interval = [tr1 tr2];
train_name = par.train_name;
tr1 = sprintf('[%.2f,%.2f]',par.train_tr1);
tr2 = sprintf('[%.2f,%.2f]',par.train_tr2);
% tr1 = par.train_tr1;
% tr2 = par.train_tr2;
% test_interval = [ts1 ts2];
test_name = par.test_name;
ts1 = sprintf('[%.2f,%.2f]',par.test_ts1);
ts2 = sprintf('[%.2f,%.2f]',par.test_ts2);
% ts1 = par.test_ts1;
% ts2 = par.test_ts2;
class = par.class;
method = par.method;
n_Features = par.n_Features;

try m = par.m;
catch
    m=0;
end
try indMi = par.indMi;
catch
    indMi = 0;
end

% Accuracy
ResultAcc.date = datetime('now');
ResultAcc.method = method;
ResultAcc.subj = subj;
ResultAcc.file = file;
ResultAcc.train_name = train_name;
ResultAcc.train_num = par.elem_train;
ResultAcc.train_start = tr1;
ResultAcc.train_stop = tr2;
ResultAcc.train_Acc = mean([Acc.acc_train]);
ResultAcc.valid_name = train_name;
ResultAcc.valid_num = par.elem_valid;
ResultAcc.valid_Acc = mean([Acc.acc_valid]);
ResultAcc.test_name = test_name;
ResultAcc.test_num = par.elem_test;
ResultAcc.test_start = ts1;
ResultAcc.test_stop = ts2;
ResultAcc.test_Acc = mean([Acc.acc_test]);
ResultAcc.m = m;
for ncl=1:length(class)
    class_name = sprintf('Class%d',ncl);
    ResultAcc.(class_name) = class{ncl};
end
ResultAcc.indMi = indMi;
ResultAcc.Filter = par.Filter;
ResultAcc.attenuation = par.attenuation;
ResultAcc.TotalFeatures = par.TotalFeatures;
ResultAcc.nFeatures= n_Features;
ResultAcc.kfold_Splitdata = par.kfold;
ResultAcc.irng = par.irng;