function   [mtl,stl,bci]=CI_compute(X,par)
% function [mtl,stl,bci]=CI_compute(X,par)
try 
    P=par.P;    % Confidence interval
catch
    % P=90;	    % Default Confidence Interval
    P=95;	    % Default Confidence Interval
end
try 
    nBootSamples=par.N; % number of bootstrap samples
catch
    nBootSamples=10000;
end
try
    OPT_CALC=par.opt;
catch
    % OPT_CALC=[1,0,0];  % BOOTSTRAP
    % OPT_CALC=[0,1,0];  % BOOTSTRAP, equivalent to previous option - slower and demonstrative
    % OPT_CALC=[0,0,1];  % mean and variance
    OPT_CALC=[0,0,1];  
end
alpha =(1-P/100);
% xval=all{n,k};
xval=X;
if OPT_CALC(1)
    %% WITH BOOTSTRAP
    try
        [bci,bmeans] = bootci(nBootSamples,{@mean,xval},'alpha',alpha);
    catch
        [bci,bmeans] = bootci(nBootSamples,{@mean,xval},'alpha',alpha,'type','per');
    end
    mtl          = mean(bmeans);
    stl          = diff(bci)/2;                % confidence semi interval
elseif OPT_CALC(2)
    %% WITH BOOTSTRAP and HAND CALCULATED standard error
    Tmax=size(xval,2);
    mtl=nan(1,Tmax);
    stl=nan(1,Tmax);
    for iT=1:Tmax
        xval_cur= xval(:,iT);
        stats   = bootstrp(nBootSamples,@(x)[mean(x) std(x)],xval_cur);
        sts     = mean(stats);
        mtl(iT) = sts(1);
        stl(iT) = sts(2);
        % get bootstrap means and std
    end
    nSamples= size(xval,1);
    ts      = tinv([alpha/2  (1-alpha/2)],nBootSamples-1);% T-Score with bootstrap
    ts_     = max(ts);                                    % by construction the interval is simmetric
    stl     = stl/sqrt(nSamples);                         % standard error
    
    stl     = ts_*stl;                                    % confidence semi interval size
    bci     = [mtl-stl; mtl+stl];                         % confidence interval
elseif OPT_CALC(3)
    %% WITHOUT BOOTSTRAP
    nSamples     = size(xval,1);
    
    mtl   = mean(xval,1,'omitnan');
    stl   = std(xval,[],1,'omitnan')/sqrt(nSamples);     % standard error
    ts    = tinv([alpha/2  (1-alpha/2)],nSamples-1);     % T-Score
    ts_   = max(ts);

    stl   = ts_*stl;                              % confidence semi interval size
    bci   = [mtl-stl; mtl+stl];                   % confidence interval
end 


