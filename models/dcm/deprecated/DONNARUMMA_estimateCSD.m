function   CSD = DONNARUMMA_estimateCSD(LFP,Hz,dt,ARp,rsfactor)
% function CSD = DONNARUMMA_estimateCSD(LFP,Hz,dt,ARp,rsfactor)
try
    ARp;
catch
    ARp      =8;
end
try
    rsfactor;
catch
    rsfactor = 0.2;
end
% and estimate spectral features under a MAR model
%--------------------------------------------------------------------------

% Bayesian Multivariate Autoregressive Modelling (MAR)
mar     = spm_mar(LFP,ARp);
% Get spectral estimates from MAR model
mar     = spm_mar_spectra(mar,Hz,1/dt);
% Get cross spectral densities from MAR model
CSD     = spm_mar2csd(mar,Hz,1/dt); 
% rescale
CSD     = (CSD/max(real(CSD(:))))*rsfactor;       % rescale

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bayesian Multivariate Autoregressive Modelling (MAR)
% FORMAT [mar,y,y_pred] = spm_mar (X,p,prior,verbose)
%
% Matrix of AR coefficients are in form
% x_t = -a_1 x_t-1 - a_2 x_t-2 + ...... - a_p x_t-p
% where a_k is a d-by-d matrix of coefficients at lag k and x_t-k's are 
% vectors of a d-variate time series.
%
% X              T-by-d matrix containing d-variate time series
% p              Order of MAR model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get spectral estimates from MAR model
% FORMAT [mar] = spm_mar_spectra (mar,freqs,ns,show)
% mar   - MAR data structure (see spm_mar.m)
% freqs - [Nf x 1] vector of frequencies to evaluate spectra at
% ns    - samples per second (default: ns = 2*freqs(end))
% show  - 1 if you wish to plot estimates (default is 0)
%
% The returned mar will have the following fields specified:
%
% .P     [Nf x d x d] Power Spectral Density matrix
% .H     [Nf x d x d] Transfer Function matrix
% .C     [Nf x d x d] Coherence matrix
% .dtf   [Nf x d x d] Kaminski's Directed Transfer Function matrix
% .pve   [Nf x d x d] Geweke's proportion of variance explained
% .gew   [Nf x d x d] Geweke's frequency domain Granger causality
% .pdc   [Nf x d x d] Baccala's Partial Directed Coherence
% .L     [Nf x d x d] Phase matrix
% .f     [Nf x 1] Frequency vector
% .ns    Sample rate
%
% dtf(f,i,j) is the DTF at frequency f from signal j to signal i
% pdc(f,i,j) is the PDC at frequency f from signal j to signal i
% pve(f,i,j) is the proportion of power in signal i at frequency f that can
%            be predicted by signal j. 
% gew(f,i,j) is the Granger casuality from signal j to signal i at frequency f.
%            gew=-log(1-pev)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get spectral estimates from MAR model
% FORMAT [csd,dtf,coh,pha] = spm_mar2csd(mar,freqs,ns)
% mar   - MAR coefficients or structure (see spm_mar.m)
% freqs - [Nf x 1] vector of frequencies to evaluate spectra at
% ns    - samples per second
%
% csd   - cross spectral density
% coh   - coherence
% pha   - phase
% mtf   - modulation transfer function
