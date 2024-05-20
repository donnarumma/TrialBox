function [n_epoch,time_intervals,bin_intervals] = computeEpochIntervals(t_total,t_epoch,overlap_percent,fsample)

% Parametri
% t_total = par.t_total; % Durata totale in secondi
% t_epoch = par.t_epoch; % Durata di ciascuna epoca in secondi
% overlap_percent = par.overlap_percent; % Percentuale di overlap tra le epoche
% fsample = par.fsample;
% Calcola il numero massimo di epoche possibile, considerando l'overlap
n_epoch_max = floor(t_total / (t_epoch * (1 - overlap_percent / 100)));

% Calcola il numero effettivo di epoche da estrarre
n_epoch = min(n_epoch_max, floor((t_total - t_epoch) / (t_epoch * (1 - overlap_percent / 100))) + 1);

% Calcola la durata effettiva di ciascuna epoca, considerando l'overlap
effective_t_epoch = t_epoch * (1 - overlap_percent / 100);

% Inizializza la matrice degli intervalli di tempo
time_intervals = zeros(n_epoch, 2);

% Genera gli intervalli di tempo per ciascuna epoca
for i = 1:n_epoch
    start_time = max(0, (i - 1) * effective_t_epoch);
    end_time = min(t_total, start_time + t_epoch);
    time_intervals(i, :) = [start_time, end_time];
end

% Genera gli intervalli di tempo per ciascuna epoca
n_bin= floor(min((time_intervals(:,2).*fsample) - (time_intervals(:,1).*fsample)));
bin_intervals = zeros(n_epoch, 2);
for i = 1:n_epoch
    start_bin = 1 + (i - 1) * (n_bin - floor(n_bin * overlap_percent/100));
    end_bin = start_bin + n_bin - 1;
    bin_intervals(i, :) = [start_bin, end_bin];
end

% % Visualizza gli intervalli
% disp(n_epoch);
% disp(time_intervals);
% disp(bin_intervals);