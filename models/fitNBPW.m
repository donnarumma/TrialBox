function [pred,out] = fitNBPW(data_train,data_test,par)

% Naive Bayesian Parzen Window classifier

execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ', mfilename); end

labs_train              = par.labs_train;
labs_test               = par.labs_test;

typeClass       = unique(labs_train);
numClass        = length(typeClass);

nt = length(labs_test);
proba = zeros(nt,1);
index = zeros(nt,1);
ind = index;
P_out = nan(nt,numClass);
for nlabs = 1:nt
    pwx = nan(numClass,1);
    for ncl =1:length(typeClass)
        nb_test = data_test{ncl,1};
        pwx(ncl) = NBPW(data_train{ncl}, labs_train, nb_test(nlabs,:), ncl);
    end
    P_out(nlabs,:) = pwx;
    [proba(nlabs), index(nlabs)] = max(pwx);
end
for ncl = 1:numClass
    ind(index==ncl) = ncl;
end
pred = ind;
out.prob = P_out;
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end


