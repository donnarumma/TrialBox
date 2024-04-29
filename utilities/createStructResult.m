function [Result,Result_class] = createStructResult(res,par)

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
train_class = nan(length(res),length(res(1).train.Accuracy_class));
test_m = nan(length(res),1);
test_class = nan(length(res),length(res(1).test.Accuracy_class));
for i=1:length(res)
    train_m(i) = res(i).train.Accuracy;
    train_class(i,:) = res(i).train.Accuracy_class;
    test_m(i) = res(i).test.Accuracy;
    test_class(i,:) = res(i).test.Accuracy_class;
end

mean_class_train = mean(train_class,1);
std_class_train = std(train_class,0,1);
mean_class_test = mean(test_class,1);
std_class_test = std(test_class,0,1);


Result.date = datetime('now');
Result_class.Method = method;
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

% Result for class
train_cl_mean = cell(length(mean_class_train),1);
train_cl_std = cell(length(mean_class_train),1);
for i=1:length(mean_class_train)
    train_cl_mean{i} = strcat(sprintf('train_cl%d',i),'_mean');
    train_cl_std{i} = strcat(sprintf('train_cl%d',i),'_std');
end
test_cl_mean = cell(length(mean_class_test),1);
test_cl_std = cell(length(mean_class_test),1);
for i=1:length(mean_class_test)
    test_cl_mean{i} = strcat(sprintf('test_cl%d',i),'_mean');
    test_cl_std{i} = strcat(sprintf('test_cl%d',i),'_std');
end

Result_class.date = datetime('now');
Result_class.Method = method;
Result_class.subj = subj;
Result_class.file = file;
Result_class.train_name = train_name;
Result_class.train_start = tr1;
Result_class.train_stop = tr2;
for itr =1:length(train_cl_mean)
    Result_class.(train_cl_mean{itr}) = mean_class_train(itr);
    Result_class.(train_cl_std{itr}) = std_class_train(itr);
end
Result_class.test_name = test_name;
Result_class.test_start = ts1;
Result_class.test_stop = ts2;
for its=1:length(test_cl_mean)
    Result_class.(test_cl_mean{its}) = mean_class_test(its);
    Result_class.(test_cl_std{its}) = std_class_test(its);
end
Result_class.m = m;
for ncl=1:length(class)
    class_name = sprintf('Class%d',ncl);
    Result_class.(class_name) = class{ncl};
end
Result_class.irng = par.irng;
