function filtered_data = FiltBankNine(data_eeg, par)
% filtered_data = FiltBankNine(data_eeg, par)

% FILTERBANK decomposes the EEG into multiple frequency pass bands using
% causal Chebyshev Type II or Buttherworth
%
%   DETAILS
%   A total of 9 band-pass filters are used with (4-8; 8-12; 12-16; ... ; 36-40 Hz).
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
f_1_low = 4;
f_1_high = 8;

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

%% teta
f_2_low = 8;
f_2_high = 12;

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
f_3_low = 12;
f_3_high = 16;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_3_low f_3_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_3_low f_3_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,3) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end

%% beta low
f_4_low = 16;
f_4_high = 20;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_4_low f_4_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_4_low f_4_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);

% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,4) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end

%% beta high
f_5_low = 20;
f_5_high = 24;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_5_low f_5_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_5_low f_5_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,5) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end

%% gamma Low
f_6_low = 24;
f_6_high = 28;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_6_low f_6_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_6_low f_6_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,6) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end

%% gamma Low
f_7_low = 28;
f_7_high = 32;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_7_low f_7_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_7_low f_7_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,7) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end
%% gamma Low
f_8_low = 32;
f_8_high = 36;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_8_low f_8_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_8_low f_8_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,8) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));
end
%% gamma Low
f_9_low = 36;
f_9_high = 40;
if strcmp(par.TypeFilter,'Chebyshev')
    % Chebyshev type II filter
    [z,p,k] = cheby2(order,attenuation,2*[f_9_low f_9_high]/fsamp);
elseif strcmp(par.TypeFilter,'Butterworth')
    % Butterworth filter
    [z,p,k] = butter(order,2*[f_9_low f_9_high]/fsamp,'bandpass');
end
[sos,g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);
% filtering each eeg from different channels/trials/runs
for j = 1:size(data_eeg,1)
    filtered_data(j,:,9) = filtfilt(hd.sosMatrix,hd.ScaleValues,data_eeg(j,:));

end

end

