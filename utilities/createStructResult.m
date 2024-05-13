function [ResultKappa,ResultAcc,ResultAcc_class] = createStructResult(res,par)

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
train_kappa_m = nan(length(res),1);
test_kappa_m = nan(length(res),1);
for i=1:length(res)
    train_kappa_m(i) = res(i).train.kappaValue;
    test_kappa_m(i) = res(i).test.kappaValue;
end

ResultKappa.date = datetime('now');
ResultKappa.method = method;
ResultKappa.subj = subj;
ResultKappa.file = file;
ResultKappa.train_name = train_name;
ResultKappa.train_start = tr1;
ResultKappa.train_stop = tr2;
ResultKappa.train_Kappa_mean = mean(train_kappa_m);
ResultKappa.train_Kappa_std = std(train_kappa_m);
ResultKappa.test_name = test_name;
ResultKappa.test_start = ts1;
ResultKappa.test_stop = ts2;
ResultKappa.test_Kappa_mean = mean(test_kappa_m);
ResultKappa.test_Kappa_std = std(test_kappa_m);
ResultKappa.m = m;
for ncl=1:length(class)
    class_name = sprintf('Class%d',ncl);
    ResultKappa.(class_name) = class{ncl};
end
ResultKappa.irng = par.irng;
ResultKappa.Filter = par.Filter;

% Accuracy
train_acc_m = nan(length(res),1);
train_acc_class = nan(length(res),length(res(1).train.Accuracy_class));
test_acc_m = nan(length(res),1);
test_acc_class = nan(length(res),length(res(1).test.Accuracy_class));
for i=1:length(res)
    train_acc_m(i) = res(i).train.Accuracy;
    train_acc_class(i,:) = res(i).train.Accuracy_class;
    test_acc_m(i) = res(i).test.Accuracy;
    test_acc_class(i,:) = res(i).test.Accuracy_class;
end

mean_class_train_acc = mean(train_acc_class,1);
std_class_train_acc = std(train_acc_class,0,1);
mean_class_test_acc = mean(test_acc_class,1);
std_class_test_acc = std(test_acc_class,0,1);


ResultAcc.date = datetime('now');
ResultAcc_class.Method = method;
ResultAcc.subj = subj;
ResultAcc.file = file;
ResultAcc.train_name = train_name;
ResultAcc.train_start = tr1;
ResultAcc.train_stop = tr2;
ResultAcc.train_Acc_mean = mean(train_acc_m);
ResultAcc.train_Acc_std = std(train_acc_m);
ResultAcc.test_name = test_name;
ResultAcc.test_start = ts1;
ResultAcc.test_stop = ts2;
ResultAcc.test_Acc_mean = mean(test_acc_m);
ResultAcc.test_Acc_std = std(test_acc_m);
ResultAcc.m = m;
for ncl=1:length(class)
    class_name = sprintf('Class%d',ncl);
    ResultAcc.(class_name) = class{ncl};
end
ResultAcc.irng = par.irng;
ResultAcc.Filter = par.Filter;

% Result for class
train_cl_mean = cell(length(mean_class_train_acc),1);
train_cl_std = cell(length(mean_class_train_acc),1);
for i=1:length(mean_class_train_acc)
    train_cl_mean{i} = strcat(sprintf('train_cl%d',i),'_mean');
    train_cl_std{i} = strcat(sprintf('train_cl%d',i),'_std');
end
test_cl_mean = cell(length(mean_class_test_acc),1);
test_cl_std = cell(length(mean_class_test_acc),1);
for i=1:length(mean_class_test_acc)
    test_cl_mean{i} = strcat(sprintf('test_cl%d',i),'_mean');
    test_cl_std{i} = strcat(sprintf('test_cl%d',i),'_std');
end

ResultAcc_class.date = datetime('now');
ResultAcc_class.Method = method;
ResultAcc_class.subj = subj;
ResultAcc_class.file = file;
ResultAcc_class.train_name = train_name;
ResultAcc_class.train_start = tr1;
ResultAcc_class.train_stop = tr2;
for itr =1:length(train_cl_mean)
    ResultAcc_class.(train_cl_mean{itr}) = mean_class_train_acc(itr);
    ResultAcc_class.(train_cl_std{itr}) = std_class_train_acc(itr);
end
ResultAcc_class.test_name = test_name;
ResultAcc_class.test_start = ts1;
ResultAcc_class.test_stop = ts2;
for its=1:length(test_cl_mean)
    ResultAcc_class.(test_cl_mean{its}) = mean_class_test_acc(its);
    ResultAcc_class.(test_cl_std{its}) = std_class_test_acc(its);
end
ResultAcc_class.m = m;
for ncl=1:length(class)
    class_name = sprintf('Class%d',ncl);
    ResultAcc_class.(class_name) = class{ncl};
end
ResultAcc_class.irng = par.irng;
ResultAcc_class.Filter = par.Filter;
