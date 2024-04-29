function filtered_data = FiltBankEEGbands(data_eeg, par)
% filtered_data = FiltBankEEGbands(data_eeg, par)

% FILTERBANK decomposes the EEG into multiple frequency pass bands using
% causal Chebyshev Type II or Buttherworth
%
%   DETAILS
%   A total of 7 band-pass filters are used (delta,teta,alpha,beta_low,beta_high,gamma_low,gamma_high)
%   with (f_min-4; 4-8; 8-13; 13-20; 20-30; 30-40; 40-f_max Hz).
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
order = par.order;                 % * chosen arbitrarily
attenuation = par.attenuation;           % * only for Chebyshevchev
fsamp = par.fsample;

% filtered signal init
filtered_data = zeros([size(data_eeg) 7]);


%% delta
f_delta_low = par.f_min;
f_delta_high = 4;

if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_delta_low f_delta_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_delta_low f_delta_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,1) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end

%% teta
f_teta_low = 4;
f_teta_high = 8;

if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_teta_low f_teta_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_teta_low f_teta_high]/fsamp,'bandpass');
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

%% beta high
f_betaH_low = 20;
f_betaH_high = 30;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_betaH_low f_betaH_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_betaH_low f_betaH_high]/fsamp,'bandpass');
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

%% gamma Low
f_gammaH_low = 40;
f_gammaH_high = par.f_max;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_gammaH_low f_gammaH_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_gammaH_low f_gammaH_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,7) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));

end

end
