function filtered_data = FiltBankVar(data_eeg, par)
% filtered_data = FiltBankVar(data_eeg, par)

% FiltBankVar decomposes the EEG into multiple frequency pass bands using
% causal Chebyshev Type II or Buttherworth Filters 
%
%   DETAILS
%   A total of n band-pass filters are used. the n number of filter are
%   related on frequency span (f_max and f_min), frequency window and delta frequency interval.
%
%   INPUT
%   'data_eeg' is a single electroencephalographic signal in time domain;
%
%   OUTPUT
%   'filtered_eeg' is an array with eeg filtered at different pass-bands.
%
%   UPDATE: 2024/04/29
%

% filter parameters
order   = par.order;                 % * chosen arbitrarily
attenuation = par.attenuation;       % * only for Chebyshevchev
fsamp   = par.fsample;
f_min   = par.f_min;
f_max   = par.f_max;
fw      = par.fw;
delta_f = par.deltaf;

num_filters = floor((f_max - f_min -fw) / delta_f) + 1;

% filtered signal init
filtered_data = zeros([size(data_eeg) num_filters]);


%% delta
f_filt_low = f_min;
f_filt_high = f_min + fw;
n_fil = 1;
while f_filt_high < f_max
    if strcmp(par.TypeFilter,'Chebyshev')
        % Chebyshev type II filter
        [z,p,k] = cheby2(order,attenuation,2*[f_filt_low f_filt_high]/fsamp);
    elseif strcmp(par.TypeFilter,'Butterworth')
        % Butterworth filter
        [z,p,k] = butter(order,2*[f_filt_low f_filt_high]/fsamp,'bandpass');
    end
    [sos,g] = zp2sos(z,p,k);
    hd = dfilt.df2sos(sos,g);
    % filtering each eeg from different channels/trials/runs
    for j = 1:size(data_eeg,1)
        filtered_data(j,:,n_fil) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
    end
    n_fil = n_fil + 1;
    f_filt_low = f_filt_low + delta_f;
    f_filt_high = f_filt_high + fw;
end

