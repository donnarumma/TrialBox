function   [data_trials, out]=faEncode(data_trials,par)
% function [data_trials, out]=faEncode(data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
nTrials     = length(data_trials);
InField     = par.InField;
OutField    = par.OutField;

% from [nChannels,nTimes]  = size(data_trials(1).(InField));
% data are put in nChannels x nTimes * nTrials 
X_data      = [data_trials.(InField)];          
X_data      = X_data - repmat(par.mu,1,size(X_data,2));

Wfa         = params.Wfa;   % nChannels x nComponents
Ph          = params.Ph;

[nChannels, nTimesxnTrials] = size(X_data);
nComponents  = size(Wfa, 2);

XcXc      = X_data * X_data';

Ident     = eye(nComponents);

const   = -nChannels/2*log(2*pi);

iPh     = diag(1./Ph);
iPhL    = iPh * Wfa;    
MM      = iPh - iPhL / (Ident + Wfa' * iPhL) * iPhL';
beta    = Wfa' * MM;            % nComponents x nChannels

Zfa     = beta * Xc;            % nComponents x nTimesxNtrials
Zcov    = Ident - beta * Wfa;   % nComponents x nComponents; same for all observations

LL      = nTimesxnTrials*const + 0.5*nTimesxnTrials*logdet(MM) - 0.5 * sum(sum(MM .* XcXc));

Zfa    = reshape(Zfa,nComponents,nTimes,nTrials);   % nComponents x nTimes x nTrials

xfld    = par.xfld;
for iTrial=1:nTrials
    data_trials(iTrial).(OutField)          = Zfa(:,:,iTrial);
    data_trials(iTrial).([xfld OutField])   = data_trials(iTrial).([xfld InField]);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Time Elapsed: %.2f s\n',out.exectime); end
out.cov             = Zcov;
out.numComponents   = nComponents;
out.likelihood      = LL;