% In generale, la struttura dei dati per soggetto deve essere con almeno i
% campi:
% - Condizione (ad ora trialTypeCond)
% - Direzione   (ad ora trialTypeDir) 
% - Spikes (in cui ogni cella è di dimensione Neuroni*#bin_trial)
%   N.B.! i trials per questo tipo di analisi devono avere tutti medesima
%   lunghezza
% - MAnifold: strutturato come i trials
%
%A_struct=K_Struct'
%B_struct=S_Struct'

% per far girare la funzione Hasson (compute_LL_matrix): 
% intanto genero la struttura di configurazione cui passo tutte le info necessarie

config_=struct;

config_.corr_obj= "manifold"  % "manifold" | "hat-obs" | "hat-hat"
    % a seconda degli oggetti che voglio correlare, scelgo
    % manifold: correli le manifold e basta
    % Negli altri due casi, si passa per un decoder per avere dei valori
    % stimati dell'attività neurale, a partire dalle manifold stesse: 
    % - "hat-obs" correlo le stime di un soggetto con i dati osservati
    %    dell'altro (e viceversa). 
    % -  Nel caso "hat-hat" correlo le stime fatte con decoder dei soggetti

% use_neural indica se i dati neurali vengono utilizzati nella pipeline
% (non necessariamente per la correlazione). 
config_.use_neural=1

% Ci deve essere coerena tra l'uso di dati neurali ed gli oggetti da
% correlare; se sceglo manifold i dati neurali non occorrono per la
% correlazione (possono cmq occorrere per altro). 
    switch config_.corr_obj
        case "manifold"
            % use_neural can be true or false (future-proof)
        case {"hat-obs", "hat-hat"}
            assert(config_.use_neural == true, ...
           'use_neural must be true for corr_obj = %s', config_.corr_obj);
        otherwise
            error("Unknown corr_obj");
    end


 % parametri decoder (nei casi hat-obs e hat-hat, devo scegliere un
 % decoder per fare le correlazioni tra dati osservati e stimati)
if ismember(config_.corr_obj, ["hat-obs","hat-hat"])

    % struttura coerente col decoder (per ora abbiamo solo ridge reg e knn)
    config_.decoder = struct;
    config_.decoder.method = "knn";   % "ridge" | "knn" | ...

    switch config_.decoder.method
        case "ridge"
            config_.decoder.ridge_lambdas = logspace(-4,4,10);
        case "knn"
            config_.decoder.knn_k = 5;
        otherwise
            error("Unknown decoder method");
    end
end


% I tre parametri seguenti determinano il numero di blocchi (lag).
% block_stride è una quota di block_size e definisce lo SHIFT del blocco:
%   stride_abs = floor(block_stride * block_size)
% L’overlap implicito è:
%   block_size - stride_abs
% Il numero di blocchi è:
%   n_blocchi = floor((trial_length - block_size) / stride_abs) + 1
% trial length (in bins)
config_.trial_length= 164 % lunghezza dei trial (int Fix)
% trial lenght timespan in ms 
config_.time_length=500
config_.block_size =  50  % dimensione dei blocchi (lag) (int < trial_length)
config_.block_stride = 0.2  

% i due parametri successivi occorrono (Francesco docet) a determinare il
% numero di trial permettendo una media degli stessi:
% sub_block_size è la finestra temporale su cui si calcola la media locale.
% sub_block_stride controlla lo shift tra sottoblocchi (overlap implicito).
% Caso degenere:
%   sub_block_size = block_size --> una sola osservazione per trial per blocco.
%   per direzione 
config_.sub_block_size =  25  % 
config_.sub_block_stride = 1  % 

% coerenza dei parametri temporali
assert(config_.block_size < config_.trial_length, ...
    'block_size must be smaller than trial_length');
assert(config_.block_stride > 0 && config_.block_stride <= 1, ...
    'block_stride must be in (0,1]');
assert(config_.sub_block_size <= config_.block_size, ...
    'sub_block_size must be <= block_size');
assert(config_.sub_block_stride > 0 && config_.sub_block_stride <= 1, ...
    'sub_block_stride must be in (0,1]');

if config_.sub_block_size == config_.block_size && config_.sub_block_stride ~= 1
    warning(['sub_block_size == block_size: sub_block_stride has no effect ', ...
             '(single average per block).']);
end

config_.max_dir=8
config_.max_cond=3

% soggetti, direzione e condizione  (task)
% soggetto attivo e passivo dipendono dalla condiz.
results={}
for i =1:config_.max_cond
    config_.cond = [i] % condizione
    cond=i
    for j=1:config_.max_dir
        config_.dir = [j]  % list of direction(s)
        dir_=j
        %%% Just notice that:
        % In compute_LL_matrix(Y_subject, X_subject, ...),
        % Y_subject maps to the ROWS (Y-axis) of the lag–lag matrix,
        % X_subject maps to the COLUMNS (X-axis) of the lag–lag matrix.
        Y_subject=S_Struct
        X_subject=K_Struct
        % define which subject 
        config_.y_subject='s' 
        config_.x_subject='k'
        L_L_matrix = compute_LL_matrix(Y_subject,X_subject, config_);
        LL_plot(L_L_matrix,config_)

        output_dir = fullfile(pwd,'results_heatmap_manifold')
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        filename = sprintf('cond_%d_dir_%d.png', config_.cond, config_.dir);
        saveas(gcf, fullfile(output_dir, filename));
        %% Calcolo medie varie
        corr_matrix_ = L_L_matrix;
        n_lags = size(corr_matrix_, 1);
        % Maschere per isolare le parti di interesse
        % Nota: dopo inversione asse Y nei grafici:
        % triu (alto dx in mat) = sotto la diagonale nel grafico = S anticipa K (S→K)
        % tril (basso sx in mat) = sopra la diagonale nel grafico = K anticipa S (K→S)
        S2K_triangle = triu(corr_matrix_, 1);   
        K2S_triangle = tril(corr_matrix_, -1);  
        diag_values  = diag(corr_matrix_);

        % clean small values
        eps = 1e-6;
        S2K_triangle(abs(S2K_triangle) < eps) = NaN;
        K2S_triangle(abs(K2S_triangle) < eps) = NaN;
        diag_values(abs(diag_values) < eps)   = NaN;

        % MEDIE GLOBALI
        mean_total   = mean(corr_matrix_(:), 'omitnan');
        mean_S2K     = mean(S2K_triangle(:), 'omitnan');   % S guida K
        sd_S2K       = std(S2K_triangle(:), 'omitnan');
        mean_K2S     = mean(K2S_triangle(:), 'omitnan');   % K guida S
        sd_K2S       = std(K2S_triangle(:), 'omitnan');
        mean_diag    = mean(diag_values, 'omitnan');       % sincronia
        sd_diag      = std(diag_values(:), 'omitnan');

        % MEDIE PER RIGA (asse Y / S)
        mean_S2K_rows = mean(S2K_triangle, 2, 'omitnan');  % per ogni lag di S
        mean_K2S_rows = mean(K2S_triangle, 2, 'omitnan');

        % MEDIE PER COLONNA (asse X / K)
        mean_S2K_cols = mean(S2K_triangle, 1, 'omitnan');  % per ogni lag di K
        mean_K2S_cols = mean(K2S_triangle, 1, 'omitnan');

        results(cond).cond(dir_).mean_total     = mean_total;
        % globali
        results(cond).cond(dir_).mean_S2K       = mean_S2K;
        results(cond).cond(dir_).sd_S2K         = sd_S2K;
        results(cond).cond(dir_).mean_K2S       = mean_K2S;
        results(cond).cond(dir_).sd_K2S         = sd_K2S;
        results(cond).cond(dir_).mean_diag      = mean_diag;
        results(cond).cond(dir_).sd_diag        = sd_diag;

        % vettori per righe (asse Y = S)
        results(cond).cond(dir_).mean_S2K_rows  = mean_S2K_rows;
        results(cond).cond(dir_).mean_K2S_rows  = mean_K2S_rows;

        % vettori per colonne (asse X = K)
        results(cond).cond(dir_).mean_S2K_cols  = mean_S2K_cols;
        results(cond).cond(dir_).mean_K2S_cols  = mean_K2S_cols;
    end
end

%% Barplots tutti i risultati BARPLOTS 
% results=load('results_01_12_25_0_400_ms.mat')
cc = config_.max_cond;
dd = config_.max_dir;
direzioni = 1:dd;

mean_upper_all=zeros(cc,dd);
mean_lower_all=zeros(cc,dd);
mean_diag_all=zeros(cc,dd)
sd_upper_all=zeros(cc,dd)
sd_lower_all=zeros(cc,dd)
sd_diag_all=zeros(cc,dd)
%%% container per differene
diffs = zeros(cc,dd);
diags = zeros(cc,dd);

for c = 1:cc
    for d = 1:dd
        elem = results(c).cond(d);
      % medie
        mean_upper_all(c,d) = elem.mean_K2S;   % K → S
        mean_lower_all(c,d) = elem.mean_S2K;   % S → K
        mean_diag_all(c,d)  = elem.mean_diag;  % S ↔ K

        % deviazioni standard
        sd_upper_all(c,d) = elem.sd_K2S;
        sd_lower_all(c,d) = elem.sd_S2K;
        sd_diag_all(c,d)  = elem.sd_diag;

        % differenze e diagonali 
        diffs(c,d) = elem.mean_K2S - elem.mean_S2K;
        diags(c,d) = elem.mean_diag;
    end
end

figure
subplot(3,1,1);
Y_u=mean_upper_all'
sd_u=sd_upper_all'
h1=bar(Y_u, 'grouped');
hold on 
% Errorbar (opzionale)
% for i = 1:numel(h1)                 
%     x = h1(i).XEndPoints;           
%     %errorbar(x, Y_u(:,i), sd_u(:,i), sd_u(:,i), ...
%      %        'k', 'linestyle','none', 'LineWidth',1, 'CapSize',6);
% end
%h2=errorbar(sd_upper_all' ,'')
title('K → S');
ylabel('Corr');
xticklabels(arrayfun(@(x) sprintf('Dir %d', x), direzioni, 'UniformOutput', false));
%legend({'Cond 1','Cond 2','Cond 3'}, 'Location','best');
%ylim([min(mean_lower_all(:)) max(mean_upper_all(:))]); grid on;
ylim([-0.4 0.4]); grid on;
%ylim([min(Y_u(:)-sd_u(:)) max(Y_u(:)+sd_u(:))] + [-0.02 0.02]); grid on;

subplot(3,1,2);
Y_l=mean_lower_all'
sd_l=sd_lower_all'
h2=bar(Y_l, 'grouped');
hold on 
% Errorbar (opzionale)
% for i = 1:numel(h2)
%     x = h2(i).XEndPoints;
%     %errorbar(x, Y_l(:,i), sd_l(:,i), sd_l(:,i), 'k','linestyle','none','LineWidth',1,'CapSize',6);
% end
title('S \rightarrow K');
ylabel('Corr');
xticklabels(arrayfun(@(x) sprintf('Dir %d', x), direzioni, 'UniformOutput', false));
%legend({'Cond 1','Cond 2','Cond 3'}, 'Location','best');
%ylim([min(mean_lower_all(:)) max(mean_upper_all(:))]); grid on;
ylim([-0.4 0.4]); grid on;

subplot(3,1,3);
Y_d=mean_diag_all'
sd_d=sd_diag_all'
h3=bar(Y_d, 'grouped');
hold on
% Errorbar (opzionale)
% for i = 1:numel(h3)
%     x = h3(i).XEndPoints;
%     %errorbar(x, Y_d(:,i), sd_d(:,i), sd_d(:,i), 'k','linestyle','none','LineWidth',1,'CapSize',6);
% end
title('S ↔ K');
ylabel('Corr');
xticklabels(arrayfun(@(x) sprintf('Dir %d', x), direzioni, 'UniformOutput', false));
%legend({'Cond 1','Cond 2','Cond 3'}, 'Location','best');
%ylim([min(mean_lower_all(:)) max(mean_upper_all(:))]); grid on;
%ylim([min(Y_d(:)-sd_d(:)) max(Y_d(:)+sd_d(:))] + [-0.02 0.02]); grid on;
ylim([-0.4 0.4]); grid on;

legend(h1, {'Act S - Obs K','Act K - Obs S','Joint Task'}, ...
    'Orientation','horizontal', ...
    'Position',[0.35, 0.01, 0.3, 0.05]);  % [x, y, larghezza, altezza]

%% Barplot differenze
diffs_plot = diffs';
epsThr = 0.05;
isNS = abs(diffs_plot) < epsThr;

signlab = strings(dd,cc);
signlab(isNS) = "n.s.";
signlab(~isNS & diffs_plot > 0) = "K → S";
signlab(~isNS & diffs_plot < 0) = "S → K";
signlab=fliplr(signlab)
condlabels = { ...
    'C1: S act - K obs', ...
    'C2: K act - S obs', ...
    'C3: S act - K act' };

% labels = { ...
%     'Cond 1: S act - K obs', ...
%     'Cond 2: K act - S obs', ...
%     'Cond 3: S act - K act' ...
% };
% 
% % Ora creiamo le label in base a n_cond
% % condlabels = labels(1:3);
% % condlabels = flip(condlabels);    

figure; clf;
tiledlayout(ceil(dd/2),2,'TileSpacing','compact','Padding','compact');

for d = 1:dd
    nexttile; hold on;

    data = diffs_plot(d,:);
    data = flip(data)
    barh(data);
    h = barh(data, 'FaceAlpha', 0.95);
    yticks(1:cc);
    yticklabels(flip(condlabels));
    xlabel('Corr Diff');
    xlim([-0.5 0.5]); grid on;
    title(sprintf('dir %d',d));

    for c = 1:cc
        %cambia carattere
        text(0.3,c,signlab(d,c),'FontWeight','bold');
    end

     % Colori: verde >0, rosso <0, grigio per n.s.
    cd = repmat([0.7 0.7 0.7], cc, 1);  % parte tutto grigio [n_cond x 3]
    idxPos = ~isNS(d,:) & diffs_plot(d,:) > 0;  % barre significative positive
    idxNeg = ~isNS(d,:) & diffs_plot(d,:) < 0;  % barre significative negative
    cd(idxPos,:) = repmat([0.2 0.6 0.2], sum(idxPos), 1);  % verde
    cd(idxNeg,:) = repmat([0.8 0.2 0.2], sum(idxNeg), 1);  % rosso
    h.FaceColor = 'flat';      
    h.CData     = cd;         
    h.CData     = flip(cd,1);

    % % Etichette a DESTRA 
    % for d = n_cond:-1:1
    %     % posizione a destra del grafico, indipendente dal valore
    %     % aggiunge uno 0.qcs limiti
    %     xpos = max(xlim) + 0.01;  
    %     %scrive lsetichetta sulla riga d, allineata in mezzo verticalmente,
    %     % e “appoggiata” a sinistra orizzontalmente (verso l%interno del
    %     % grafico).
    %     text(xpos, d, signlab(d,c), ...
    %     'VerticalAlignment', 'middle', ...
    %     'HorizontalAlignment', 'left', ...
    %     'FontWeight', 'bold');
    % end

    hold off;
end

diags_plot = diags';

figure;
tiledlayout(ceil(dd/2),2,'TileSpacing','compact','Padding','compact');

for d = 1:dd
    nexttile; hold on
    data = diags_plot(d, :);
    data = flip(data);

    % barra
    h = barh(data, 'FaceAlpha', 0.95);

    title(sprintf('dir %d', c));
    yticks(1:cc);
    yticklabels(flip(condlabels));
    xlim([-0.5 0.5]); grid on;
    title(sprintf('dir %d',d));
end
