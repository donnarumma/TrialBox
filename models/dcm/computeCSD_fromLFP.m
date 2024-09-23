function   CSD=computeCSD_fromLFP(LFP,Hz,dt,ARp,rsfactor)
% function Trials=computeCSD_fromLFP(Trials)
% t=tic;
% ARp          = 8; 
% rsfactor     = 0.2;
% Hz           = 1:64;
NTrials=length(LFP);
CSD          = cell(1,NTrials);

for indTrial     = 1:NTrials
    % fprintf('Selecting Trial %g\n',indTrial);
    CSD{indTrial} = DONNARUMMA_estimateCSD(LFP{indTrial},Hz,dt,ARp,rsfactor);
end
% fprintf('Computed Cross Spectral Densities in %g s\n',toc(t));
