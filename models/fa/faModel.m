function   [data_trials, out]=faModel(data_trials,par)
% function [data_trials, out]=faModel(data_trials,par)

execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

tol         = par.tol;
nIters      = par.nIters;
InField     = par.InField;
nComponents = par.numComponents;
typ         = par.type;
minVarFrac  = par.minVarFrac;
% original data are nChannels x nTimes -> transformed in nChannels x nTimes*nTrials 
X_data      = [data_trials.(InField)]; 
% Wfa: nChannels x K components. 
% Zfa: nTimes*nTrials x K components
mu                         = mean(X_data,2); % nChannels x 1

[nChannels, nTimesxTrials] = size(X_data);

% Initialization of parameters
covX                = cov(X_data', 1);
if rank(covX) == nChannels
    scale           = exp(2*sum(log(diag(chol(covX))))/nChannels);
else
    % covX may not be full rank because nTimesxTrials < nChannels
    % fprintf('WARNING in fastfa.m: Data matrix is not full rank.\n');
    rankX           = rank(covX);
    eigsSorted      = sort(eig(covX), 'descend');
    scale           = geomean(eigsSorted(1:rankX));
end
Wfa         = randn(nChannels,nComponents)*sqrt(scale/nComponents);
Ph          = diag(covX);
varFloor    = minVarFrac * diag(covX);  

Ident       = eye(nComponents);
const       = -nChannels/2*log(2*pi);
LLi         = 0; 
LL          = [];
    
for iter = 1:nIters
    % =======
    % E-step
    % =======
    iPh         = diag(1./Ph);
    iPhL        = iPh * Wfa;  
    MM          = iPh - iPhL / (Ident + Wfa' * iPhL) * iPhL';
    beta        = Wfa' * MM;    % nComponents x nChannels
    
    covX_beta   = covX * beta'; %   nChannels x nComponents
    EZZ         = Ident - beta * Wfa + beta * covX_beta;
        
    % Compute log likelihood
    LLold   = LLi;    
    ldM     = sum(log(diag(chol(MM))));
    LLi     = nTimesxTrials*const + nTimesxTrials*ldM - 0.5*nTimesxTrials*sum(sum(MM .* covX)); 
    if verbose
        fprintf('EM iteration %5i lik %8.1f \r', iter, LLi);
    end
    LL      = [LL LLi];    
    
    % =======
    % M-step
    % =======
    Wfa     = covX_beta / EZZ;
    Ph      = diag(covX) - sum(covX_beta .* Wfa, 2);
    
    if isequal(typ, 'ppca')
        Ph  = mean(Ph) * ones(nChannels, 1);
    end
    if isequal(typ, 'fa')
        % Set minimum private variance
        Ph  = max(varFloor, Ph);
    end
       
    if iter<=2
        LLbase = LLi;
    elseif (LLi < LLold)
        disp('VIOLATION');
    elseif ((LLi-LLbase) < (1+tol)*(LLold-LLbase))
        break;
    end
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
out.Wfa           = Wfa;
out.mu            = mu;
out.numComponents = numComponents;
out.likelihood    = LL;