function   [data_trials,out]=csdCompute(data_trials,par)
% function [data_trials,out]=csdCompute(data_trials,par)

execinfo=par.exec;
if ~isempty(execinfo); t=tic; end
InField     = par.InField;
OutField    = par.OutField;
ARp         = par.ARp;      % Order of MAR model
rsfactor    = par.rsfactor; % rescale factor
Hz          = par.Hz;
optrescale  = par.optrescale; % to be removed;
xfld        = 't';
NTrials     = length(data_trials);
% rsfactor    = 0.2;
% dt           = mean(diff(data_trials(1).([xfld InField]))); % sample time
dt          = 1/(2*Hz(end));
for itrial  = 1:NTrials
    % fprintf('Selecting Trial %g\n',itrial);
    timein  = data_trials(itrial).([xfld InField]);
    %dt      = mean(diff(timein));
    LFP     = data_trials(itrial).(InField)'; % time length x number of sources
    % CSD = DONNARUMMA_estimateCSD(LFP,Hz,dt,ARp,rsfactor);
    % Bayesian Multivariate Autoregressive Modelling (MAR)
    mar                                     = spm_mar(LFP,ARp);
    % Get spectral estimates from MAR model
    mar                                     = spm_mar_spectra(mar,Hz,1/dt);
    % Get cross spectral densities from MAR model
    CSD                                     = spm_mar2csd(mar,Hz,1/dt); % nFreq x nSources x nChannels
    % RESCALE
    switch optrescale
        case 1
            % divide per max and rescale 20 to max
            rescale                         = max(real(CSD(:)))/rsfactor;       % rescale
        case 2
            rescale                         = 5*max(real(CSD(:)));
        case 3 % spm_csd
            rescale                         = sqrt(sum(abs(CSD(:)).^2));
        case 4
            rescale                         = 4*std(CSD(:)); 
        case 5
            rescale                         = (max(abs(CSD(:)) - min(abs(CSD(:)))))/4;
        case 6 
            rescale                         = (max(abs(CSD(:)) - min(abs(CSD(:)))));
    end
    
    if optrescale==-1
        data_trials(itrial).(OutField)      = (CSD/max(real(CSD(:))))*rsfactor;       % rescale original
    end

    if optrescale>0
        CSD                                 = CSD/rescale;
    end
    data_trials(itrial).(OutField)          = CSD; %CSD/rescale;
    
    %fprintf('Warning\n');
    
    data_trials(itrial).([xfld OutField])   = Hz;
end


if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end