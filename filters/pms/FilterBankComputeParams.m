function   par = FilterBankComputeParams(par)
% function par = FilterBankComputeParams(par)

par.order  = 4; % Filter Order
par.attenuation = 40; % Filter Attenuation (only if it is applicable)
par.f_min  = 0.5; % min frequency range in Hz
par.f_max  = 40.5; % Max frequency range in Hz
par.fw     = 4; % frequency interval in HZ (only for FilTBankVar)
par.deltaf = 2; % sliding window frequency interval in HZ (only for FilTBankVar)
par.TypeFilter = 'Chebyshev'; % user can choose: 'Chebishev' or 'Butterworth'
par.FilterBank = 'One'; % user can choose: 
                               % 'One' (FiltBankOne.m) 
                               % 'EEGbands' (FiltEEGbands.m)
                               % 'Nine' (FiltNine.m) see mat file
                               % or 'Var' (FilTBankVar.m) see mat file for
                               % definition
par.InField         = 'y';
par.OutField        = 'y';
    try
        fnames=fieldnames(parSource);
        for in=1:length(fnames)
            par.(fnames{in})=parSource.(fnames{in});
        end
    catch
    end
end