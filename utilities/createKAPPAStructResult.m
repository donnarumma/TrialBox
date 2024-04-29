function Result = createKAPPAStructResult(res,par)

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
train_m = nan(length(res),1);
test_m = nan(length(res),1);
for i=1:length(res)
    train_m(i) = res(i).train.kappaValue;
    test_m(i) = res(i).test.kappaValue;
end

Result.date = datetime('now');
Result.method = method;
Result.subj = subj;
Result.file = file;
Result.train_name = train_name;
Result.train_start = tr1;
Result.train_stop = tr2;
Result.train_Acc_mean = mean(train_m);
Result.train_Acc_std = std(train_m);
Result.test_name = test_name;
Result.test_start = ts1;
Result.test_stop = ts2;
Result.test_Acc_mean = mean(test_m);
Result.test_Acc_std = std(test_m);
Result.m = m;
for ncl=1:length(class)
    class_name = sprintf('Class%d',ncl);
    Result.(class_name) = class{ncl};
end
Result.irng = par.irng;