function filtered_data = FiltBankPrior(data_eeg, par)
% filtered_data = FiltBankPrior(data_eeg, par)

% FILTERBANK decomposes the EEG into multiple frequency pass bands using
% causal Chebyshev Type II or Buttherworth
%
%   DETAILS A total of 6 band-pass filters are used with Alpha+Beta [8-30]
%   Hz','Gamma [40-90] Hz','Alpha [8-13] Hz','Beta Low [13-20] Hz',
%   'Beta [13-30] Hz','Gamma Low [30-40] Hz INPUT 'data_eeg' is a single
%   electroencephalographic signal in time domain;
%
%   OUTPUT
%   'filtered_eeg' is an array with eeg filtered at different pass-bands.
%
%   UPDATE: 2024/11/08
%

% filter parameters
order = par.order;                 % * chosen arbitrarily
attenuation = par.attenuation;           % * only for Chebyshevchev
fsamp = par.fsample;

% filtered signal init
filtered_data = zeros([size(data_eeg) 6]);


%% Alpha + Beta
f_1_low = 8;
f_1_high = 30;

if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_1_low f_1_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_1_low f_1_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,1) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end

%% Gamma
f_2_low = 40;
f_2_high = 90;

if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_2_low f_2_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_2_low f_2_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,2) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end

%% alpha
f_alpha_low = 8;
f_alpha_high = 13;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_alpha_low f_alpha_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_alpha_low f_alpha_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,3) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end

%% beta low
f_betaL_low = 13;
f_betaL_high = 20;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_betaL_low f_betaL_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_betaL_low f_betaL_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);

% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,4) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end

%% beta
f_gammaL_low = 13;
f_gammaL_high = 30;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_gammaL_low f_gammaL_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_gammaL_low f_gammaL_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,5) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end

%% gamma Low
f_gammaL_low = 30;
f_gammaL_high = 40;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_gammaL_low f_gammaL_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_gammaL_low f_gammaL_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,6) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end
end

