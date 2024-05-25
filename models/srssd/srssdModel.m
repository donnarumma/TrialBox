function   [data_trials, out]=srssdModel(data_trials,par)
% function [data_trials, out]=srssdModel(data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

InField     = par.InField;
% cat(3,dataTrials.(fn)) put trials on the third dimension,
X3d         = cat(3,data_trials.(InField));                     % nChannels x nTime x nTrials.
X3d         = permute(X3d,[3,1,2]);                             % nTrials x nChannels x nTimes
% each row is a "synergy" -> [x_1(1), ... x_nChannels(1), x_1(2), ... x_nChannels(2), ... , x_1(nTimes), ... x_nChannels(nTimes),    
Xdata       = reshape(X3d,size(X3d,1),size(X3d,2)*size(X3d,3)); % nTrials x nChannels*nTimes

[Manifold, Dsrssd, Wsrssd] = sr_ssd(Xdata, par.spG, par);

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
out.Manifold= Manifold; % nTrials x nAtoms
out.Wsrssd  = Wsrssd;   % nAtoms x nChannels*nTimes
% out.Dsrssd  = Dsrssd; % nChannels*nTimes x nAtoms
out.Dsrssd  = Dsrssd';  % nAtoms x nChannels*nTimes
