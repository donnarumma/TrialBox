function filtered_data= FiltBankOne(data_eeg, par)
% function filtered_data= FiltBankOne(data_eeg, par)

% FiltBankOne decomposes the EEG into single frequency pass bands using
% causal Chebyshev Type II filter or Butterworth.
%
%   DETAILS
%   A total of 1 band-pass filters are used (min_freq-max_freq Hz).
%
%   INPUT
%   'data_eeg' is a single electroencephalographic signal in time domain;
%
%   OUTPUT
%   'filtered_eeg' is an array with eeg filtered at different pass-bands.
%
%   UPDATE: 2024/04/29
%   filter_type select witch filter do you want to use:

% filter parameters
order = par.order;                 % * chosen arbitrarily
attenuation = par.attenuation;     % * only for Chebychev
fsamp = par.fsample;

% filtered signal init
filtered_data = zeros([size(data_eeg) 1]);

f_low = par.f_min;
f_high = par.f_max;

if strcmp(par.TypeFilter,'Chebyshev')

    [z,p,k] = cheby2(order,attenuation,2*[f_low f_high]/fsamp);

elseif strcmp(par.TypeFilter,'Butterworth')

    [z,p,k] = butter(order,2*[f_low f_high]/fsamp,'bandpass');
end

[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);

% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,1) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end
end