%% Ref libreria python https://zenodo.org/records/12359299
% parametrizzare blocchi
% label temporale
% distinzione 
% Replica moooolto semplificata del paper
% mi sono focalizzato sui punti chiave; 
% 1) usiamo i dati di s e k nelle reciproche condizioni act e obs;
% 2) I dati sono ricampionati (binnati)...nel senso che ogni trial,
% originarimente di lunghezza M ms, viene ridotto ad una frazione dello
% stesso (a seconda di come viene la manifold). Il ricampionamento avviene
% sommando i dati in mini finestre e mediando per la dimensione della
% finestra     
% 3)  Genero gli embeddidng con CEBRA behaviour (diciamo supervised)
%     usando più o meno sempre la stessa struttura: 
%  
% 4) ho allineato i dati campionando come fanno Husson e Zada
% 3) regressione ridge per lag e per elettrodo (con opione per  cross 
% validation esterna e nested per ricerca lambda ottimale)
% 5) calcolo correlazioni tra valori predetti di una scimmia e osservati
% dell'altra e viceversa

% tutto è più o meno parametrizzato, quello che va specificato lungo il
% codice sono:
% la condizione di interesse e quale soggetto sia arrivo o passivon nella
% specifica condizione
config.condition =  1;
config.active_subject = 's';
config.obs_subject = 'k';

% (i dati di ricampionamento usati nell'encoder)
% ampiezza delle finestre
window_=10 ;
% stride (overlap=window-shift_)
shift_=2 ;

% trials per tutte le direzioni
n_trials=64; 
% lunghezza del (singolo) trial originaria 
trial_length=400;
% la direzione che ci interessa specificamente:  filtro per direzione (1-8)
dir_=8
% il numero di trials per la direzione di interesse
n_trial=8

% la dimensione dei blocchi (lag all'interno dei trial e lo shift per
% ovrlap)
block_size=50
shift_block=20

% Nuovi PArametri: sotto-medie  per ogni blocco...invece che mediare su
% tutto il  blocco, mediamo su più punti all'intenro dello stesso creando di 
% di fatto un dato longitudinale (ho più punti all'interno di ogni trial)
% rolling windows di dimensione sub_Block, sub_Stride
sub_block=30
% sub_Stride=sub_block --> no overlap
sub_stride=10
% caso degenere:  sub_block=length block (block size)-sub_Stride=0


%% creo struttura per accesso dinamico ai dati 
% il file cebra results contiene i dati neurali, gli embedding, il
% conteggio dei trial e le label (direzioni) corrispettive. Ovviamente
% tutto è resamplato (cfr.  
load("cebra_results.mat")

DATA.s.cond_1.neural_active = s_cond_1_neural_active;
DATA.s.cond_1.label_active = s_cond_1_label_active';
DATA.s.cond_1.embed_active = s_cond_1_embed_active;
DATA.k.cond_1.neural_obs = k_cond_1_neural_passive;
DATA.k.cond_1.label_obs = k_cond_1_label_passive';
DATA.k.cond_1.embed_obs = k_cond_1_embed_passive;

DATA.k.cond_2.neural_active = k_cond_2_neural_active;
DATA.k.cond_2.label_active = k_cond_2_label_active';
DATA.k.cond_2.embed_active = k_cond_2_embed_active;
DATA.s.cond_2.neural_obs = s_cond_2_neural_passive;
DATA.s.cond_2.label_obs = s_cond_2_label_passive';
DATA.s.cond_2.embed_obs = s_cond_2_embed_passive;

DATA.k.cond_3.neural = k_cond_3_neural_joint;
DATA.k.cond_3.label = k_cond_3_label_joint';
DATA.k.cond_3.embed = k_cond_3_embed_joint;
DATA.s.cond_3.neural = s_cond_3_neural_joint;
DATA.s.cond_3.label = s_cond_3_label_joint';
DATA.s.cond_3.embed = s_cond_3_embed_joint;


%%  DATI Neuralii ed embedding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTE SONO LE INFO INIZIALI (Mutevoli)
% numero trial e conseguente  Numero di labels 

% 64 trials (totali, 8 per direzione)  da T (variabile) ms 
% (e una osservazione per ms)

% lunghezza trial resamplati
trial_len_res=floor((trial_length-window_)/shift_);

% durata in secondi (64 trials da 400 ms...ridotti a  bin)
% circa 25" totali (se considero tutti i trial nella loro interezza)
dd_total= n_trials*trial_length/1000; 
% frequenza campionamento (numero di punti al secondo)
fs=ceil(1000/(trial_length/trial_len_res));

%% Carichiamo dati relativi a una specifica condizione
%(una scimmia agisce. l'altra osserva)
% condizione 1: S agisce K osserva
% Condizione 2: K agisce S osserva
% Condizione 3: S agisce K agisce
% condizioni > 3 varianti eventuali: 
% cond 4: sono tutte e 3 le condzioni precedenti messe insieme (da inserire)
%%%%%%%%%%%%%%%%%%%%%%%%% LOAD DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[active_neural, active_label, active_embed, ...
 obs_neural, obs_label, obs_embed] = get_data_from_config(config, DATA);

%%%%%%%%%%%%%%%%%%%%%%%%%% filtra per direzione %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  lunghezza trial in tempo
trial_length_=trial_length
% filtrare poi i dati partedno dalle label
label_indices=find(active_label==dir_)
% tempo per 8 trials per direzione
%  un singolo trial
tt_single = (0:trial_len_res-1) * shift_ + window_;  
% durata complessiva in bin di tutti i trials
tt = zeros(1, n_trial * trial_len_res);
for i = 1:n_trial
    start_idx = (i-1) * trial_len_res + 1;
    stop_idx  = i * trial_len_res;

    trial_start = (i-1) * trial_length_;  % in ms
    % ogni bin è una "media di momenti"
    tt(start_idx:stop_idx) = tt_single + trial_start;
end     
% tempo in secondi della somma dei trials (sempre 8 trials per driezione)
% tt secondi
tt = tt / 1000; %
% Dati neruali di partenza
% % act: la scimmia  che muove il braccio
% % obs:e' la scimmia che osserva
% Estraggo dati relativi alla direzione indicata
active_neural=active_neural(label_indices,:);
obs_neural=obs_neural(label_indices,:);
active_label=active_label(label_indices,:);
obs_label=obs_label(label_indices,:);
% manifold generate con CEBRA (o altro)
% idem come sopra
active_embed = active_embed(label_indices,:);
obs_embed = obs_embed(label_indices,:);

% punto di onset iniziale (per trial)
onset_bin=find(tt>=0,1)
%%% sostanzialmente è dove principia ongi trial
mov_onsets=tt(onset_bin:trial_len_res:end)


%% ALLINEAMENTO DATI NEURALI - DATI movimento Embedding %%%%%%%%%%%%%%%%
% % di fatto si analizza l'attività neurale in intervalli di tempo specifici
% % successivi al principio del movimento (Hasson&Co. lo fanno intorno
% % agli onsets delle parole - 4 secondi prima e 4 dopo ogni parola).
% % come loro, che creano finestre di 250ms con overlap di 62.5 ms, creiamo
% % finestre con overlap; in particolare,  considero finestre di block_size ms 
% %  con overlap di block_size-shift ms (scelte all'innizio)
% % Analogamente procediamo per manifold generate con CEBRA.
% % la scelta delle dimmensioni dei blocchi (lags) è arbitraria...
% %  cerchiamo di includere tutti i dati. 

for i = 1:n_trial
    % limiti assoluti del trial i nel vettore concatenato ---
    trial_start_idx = (i-1)*trial_len_res + 1;
    trial_end_idx   =  i   *trial_len_res;

    %  onset assoluto del trial i 
    % onset_bin è la posizione dell'onset ALL'INTERNO di UN trial
    onset_pos = trial_start_idx + (onset_bin - 1);

    %  genera blocchi da onset, passo shift_block, ampiezza block_size ---
    start_idx   = onset_pos;
    blocchi_idx = [];
    while start_idx + block_size - 1 <= trial_end_idx
        end_idx     = start_idx + block_size - 1;
        blocchi_idx = [blocchi_idx; start_idx, end_idx]
        start_idx   = start_idx + shift_block
    end

    % estrazione dati per ciascun blocco 
    nB = size(blocchi_idx, 1);
    dati_obs{i}           = cell(nB, 1);
    dati_active{i}        = cell(nB, 1);
    manif_active_match{i} = cell(nB, 1);
    manif_obs_match{i}    = cell(nB, 1);

    for j = 1:nB
        start_blocco = blocchi_idx(j, 1);
        end_blocco   = blocchi_idx(j, 2);
        dati_active{i}{j}        = active_neural(start_blocco:end_blocco, :);
        dati_obs{i}{j}           = obs_neural(start_blocco:end_blocco, :);
        manif_active_match{i}{j} = active_embed(start_blocco:end_blocco, :);
        manif_obs_match{i}{j}    = obs_embed(start_blocco:end_blocco, :);
    end
end

%%% come posso ricampionare??? 
%%%%%%%%%%%%%%%%%%%%%%%%%% RIDGE REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creo strutture dati idoneee
% L'attività neurale (Y) ha struttura cella per movimento (trial);
% ogni cella contiene blocchi/lag di   dimensione  [t_lag × #neuroni]
% N.B.! T_lag=block_size!
% Tuttavia, anche gli embedding (X) hanno dimensione temporale: ogni X è
% una cella per movimento, contenente [t_lag × dimensione_embedding]
% dovrò fae media su block_size (durata del lag).
% di fatto, quidni, mi ritrovo con una Y per ogni combinazione lag/blocco
% di dimensione n_trial*1 (n_trial=#celle per struttura per scimmia, dati_obs, 
% dati_active)


% 1) Costruisco le y e recupero le X 
% Numero di blocchi (lags) per movimento (trial-direzione)
n_blocchi = size(blocchi_idx, 1); 
% canali...per ora li tengo distinti (ognuno con la sua numerosità)
% (nel paper i canali sono gli stessi per speaker e listener)
n_channels_active=size(active_neural,2);
n_channels_obs=size(obs_neural,2);
n_output=3;
% UNICO VALORE 
% canali (unica numerosità per entrambi = min(Ch(s.k))
%chns=min(n_channels_s,n_channels_k)

% struttura Y (S e K)
% (saranno tanti vettori di dimensione (numero trials*1) 
% un vettore per ogni lag per ogni elettrodo...dal momento che  ogni lag (per
% ogni elettrodo) include un blocco di osservazioni per movimento, si fa una
% media sul blocco...quidni un valore per ogni parola per ogni lag per ogni elettrodo
% Blocchi di variabili di rispsota per soggetto

%%%%%%%%%%%%%%%%%%% SCIMMIA active
Y_active= cell(n_blocchi, n_channels_active);
X_active= cell(n_blocchi,1);
% Y_active 
% Ogni elemento Y_active{b,e} è un vettore lungo di dimensione [t_lag * n_trial × 1]
% ottenuto concatenando nel tempo le osservazioni neurali relative al neurone e
% per ciascun movimento p, nel blocco/lag b.

%%% se voglio i blocchi senza annullare la dimensione temporael devo
%%% togliere il mean all'interno del ciclo quando formo le Y_active
% ad ogni modo, quando faccio media sulla lunghezza del blocco mi ritrvo
% con 64 osservazioni (per blocco per elettrodo) che sono le medie dei
% blocchi per ogni movimento . 
% se d'altra parte non faccio la media, per ogni trial (64) per ogni blocco
% per ogni elettrodo, mi trovo con una colonna di osservazioni che ha come
% ulteriore info la lungheza del blocco che rende il vettore lungo 64
% osservaioni ognuna mpltoiplicata per la lunghezza del blocco (estesa)
% quidni se i blocchi sono da 20 obs il vettore diventa
% (n_trial*len_blocco)*1 (è una struttura longitudinale)

% ancora, posso fare delle "sottomedie" per trial (per blocco)


for b = 1:n_blocchi
    for e = 1:n_channels_active
        Y_active{b,e} = [];

        for p = 1:n_trial
            current_block = dati_active{p}{b}(:, e);
            L = length(current_block);
            for k=1:sub_stride:(L-sub_block+1)
                idx_start=k
                idx_end=idx_start+sub_block-1;
                           
                Y_active{b,e} = [Y_active{b,e}; mean(current_block(idx_start:idx_end))];
            end
        end
    end
end

% sugli embedding X
for b = 1:n_blocchi
    X_active{b} = [];

    for p = 1:n_trial
        emb_block = manif_active_match{p}{b};   
        L = size(emb_block,1);                                  

        for k=1:sub_stride:(L-sub_block+1)
                idx_start=k
                idx_end=idx_start+sub_block-1;             
                X_active{b} = [X_active{b}; mean(emb_block(idx_start:idx_end,:),1)];
            end
    end
end



%%%%%%%%%%%%%%% SCIMMIA passiva (Obs)
% stessi sub_block e sub_stride ovviamenrte
Y_obs = cell(n_blocchi, n_channels_obs);
for b = 1:n_blocchi
    for e = 1:n_channels_obs
        Y_obs{b,e} = [];

        for p = 1:n_trial
            current_block = dati_obs{p}{b}(:, e);
            L = length(current_block);
             for k=1:sub_stride:(L-sub_block+1)
                idx_start=k
                idx_end=idx_start+sub_block-1;            
                Y_obs{b,e} = [Y_obs{b,e}; mean(current_block(idx_start:idx_end))];
            end
        end
    end
end

% X_obs
X_obs = cell(n_blocchi,1);

for b = 1:n_blocchi
    X_obs{b} = [];
    for p = 1:n_trial
        emb_block = manif_obs_match{p}{b};
        L = size(emb_block,1);
        mini = floor(L / n_sub);

        for k=1:sub_stride:(L-sub_block+1)
                idx_start=k
                idx_end=idx_start+sub_block-1;
                X_obs{b} = [X_obs{b}; mean(emb_block(idx_start:idx_end,:), 1)];
        end

    end
end

%% (Ridge) Regression and correlation
% Questo di fatto è il punto 1 nel paragrafo Model-based brain-to-brain
% coupling  Metodi (pag e4)
% Stimiamo una regressione per lag e per neurone come nel paper originale
% la X contiene le manifold allineate nel tempo  per blocco temporale. 
% Dopo l'addestramento Le risposte neurali predette (y_hat) per un soggetto
% vengono correlate con le risposte osservate dell'altro. 
% RIDGE REGRESSION per stimare i parametri utili a B2B coupling 
% NO Cross validation (solo GCV per la scelta del lambda ottimo)
% ref per GCV: "The Elements of Statistical Learning " (Hatie et, al 2008) 
%          Generalized Cross-Validation as a Method for Choosing a Good
%          Ridge Parameter, (Golub et al., 1979 Technometrics)

% lag length
t_lag=block_size
n_movements=n_trial;

% Inizializzazione delle strutture per i coefficienti, predizioni e
% correlazioni  
% Scimmia active
coefs_active = cell(n_blocchi, n_channels_active);
Y_hat_active= cell(n_blocchi, n_channels_active);

best_lambda_act = zeros(n_blocchi, n_channels_active);
lambdas = logspace(-4, 4, 10);
n_lambda = numel(lambdas);

for b = 1:n_blocchi
    X_ = X_active{b};
    d = size(X_, 2);
    I = eye(d);

    for e = 1:n_channels_active
        y_ = Y_active{b, e};
        l_y = length(y_);

        % STEP 1: Selezione lambda via GCV 
        gcv_scores = zeros(n_lambda, 1);
        for i = 1:n_lambda
                lambda = lambdas(i);
                beta = (X_' * X_ + lambda * I) \ (X_' * y_);
                y_hat_ = X_ * beta;
                H_diag = sum((X_ / (X_' * X_ + lambda * I)) .* X_, 2);
                gcv_scores(i) = mean(((y_ - y_hat_) ./ (1 - H_diag)).^2);
        end
        
        [~, best_i] = min(gcv_scores);
        lambda_best = lambdas(best_i);
        best_lambda_act(b, e) = best_i;

        % Regressione finale
        beta = (X_' * X_ + lambda_best * I) \ (X_' * y_);
        coefs_active{b, e} = beta;
        Y_hat_active{b,e}= X_ * beta;
        %%% da creare contenitore prima
        corrs_active(b, e) = compute_metric(y_, X_ * beta, y_);
           
    end 
end

% come loro prendo una soglia (che loro caclolano empiricamente con
%permutazion e tengo solo le correlazioni oltre quella soglia 
% r_threshold = 0.05;  
% corr_active_ = abs(corrs_active) > r_threshold;
% corr_active_ = corrs_active .* corr_active_ ;

%%%%%%%%%%%%  Scimmia obs
coefs_obs =cell(n_blocchi, n_channels_obs);
Y_hat_obs= cell(n_blocchi, n_channels_obs);
% store predicton... poi
%y_hat_obs_cv =cell(n_blocchi, n_channels_obs);
corrs_obs=  zeros(n_blocchi, n_channels_obs);
best_lambda_obs = zeros(n_blocchi, n_channels_obs);

for b = 1:n_blocchi
    X_ = X_obs{b};
    d = size(X_, 2);
    I = eye(d);

    for e = 1:n_channels_obs
        y_ = Y_obs{b, e};
        l_y = length(y_);

               % STEP 1: Selezione lambda via GCV 
        gcv_scores = zeros(n_lambda, 1);
        for i = 1:n_lambda
                lambda = lambdas(i);
                beta = (X_' * X_ + lambda * I) \ (X_' * y_);
                y_hat_ = X_ * beta;
                H_diag = sum((X_ / (X_' * X_ + lambda * I)) .* X_, 2);
                gcv_scores(i) = mean(((y_ - y_hat_) ./ (1 - H_diag)).^2);
        end
        
        [~, best_i] = min(gcv_scores);
        lambda_best = lambdas(best_i);
        best_lambda_obs(b, e) = best_i;


        % Regressione finale
        beta = (X_' * X_ + lambda_best * I) \ (X_' * y_);
        coefs_obs{b, e} = beta;
        Y_hat_obs{b,e} = X_ * beta;
        %%% da creare contenitore prima
        corrs_obs(b, e) = compute_metric(y_, X_*beta, y_);           
    end 
end

% % come loro prendo una soglia (che loro caclolano empiricamente con permutazioni
% % e tengo solo le correlazioni oltre quella soglia--
% r_threshold = 0.05;  
% corr_obs_=  abs(corrs_obs) > r_threshold;
% corr_obs_ = corrs_obs .* corr_obs_ ;

%%
%%%%%%%%%%%%%%%%%%%%%%% 4) Brain  to Brain coupling %%%%%%%%%%%%%%%%
% punto 2 nel paragrafo Model-based brain-to-brain coupling Metodi (pag e4)
%%%%%%%%%%%%%%%%%%  GENERALIZATION ACROSS SUBJECTS ######################
% Uso  le previsioni fatte per un soggetto (y_hat) per predire l'attività 
% neurale dell'altro (y osservata)
% Questo permette di testare se il pattern neurale è condiviso tra i soggetti
% (ad es. act ricostruisce i pattern neurali di obs e viceversa).

% % dapprima devo mediare sui neuroni (canali/elettordi)
% da un lato ho i dai neurali osservati dei due soggetti
% dall'altro le y_hat (stimate/previste dal modello)
% creo quindi i dati per le correlazioni
% mediamo sugli elettrodi
Y_act_3D = cell2mat(permute(Y_active, [3, 2, 1]));   
Y_obs_3D = cell2mat(permute(Y_obs, [3, 2, 1]));   
Y_active_mean_E= squeeze(mean(Y_act_3D, 2))';
Y_obs_mean_E= squeeze(mean(Y_obs_3D, 2))';

Y_hat_act_3D = cell2mat(permute(Y_hat_active, [3, 2, 1]));   
Y_hat_obs_3D = cell2mat(permute(Y_hat_obs, [3, 2, 1]));
Y_hat_active_mean=squeeze(mean(Y_hat_act_3D,2))'
Y_hat_obs_mean=squeeze(mean(Y_hat_obs_3D,2))'

% ottengo quindi due matrici per soggetto (hat e obs) che dovranno essere 
% reciproicamente correlate hat_active con actual_obs e hat_obs con
% actual_active

%% COSTRUZIONE HEATMAPS lag-lag fig 3a e S2 %%%%%%%%%%
% Ref metodi e4, par. "Model-based brain-to-brain coupling", Caption Fig. S2

n_lags=n_blocchi;

% Soglia minima per eliminare zeri o valori troppo piccoli
epsilon = 1e-4;
% prima di tutto le matrici lag-lag che si ottengono correlando ogni riga
% delle matrici di y previste
corr_matrix_act_obs= corr(Y_hat_active_mean',Y_obs_mean_E')
% qui traspongo per avere sempre listener sulle righe e speaker sulle
% colonne
corr_matrix_obs_act= corr(Y_hat_obs_mean',Y_active_mean_E')'

%%% Per ottenere una sola matrice, loro mediano sulle 2 ottenute
corr_avg_ = (corr_matrix_act_obs + corr_matrix_obs_act') / 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% traspsizione per condizione 2
if config.condition ==2
    corr_avg=corr_avg_'
    y_lab=['Lag (ms) ',upper(config.obs_subject)];
    x_lab=['Lag (ms) ', upper(config.active_subject)];
else
    corr_avg=corr_avg_
    y_lab=['Lag (ms) ',upper(config.active_subject)];
    x_lab=['Lag (ms) ',upper(config.obs_subject)];
end


max_corr=max(corr_avg,[],"all")
min_corr=min(corr_avg,[], "all")
%tick_indices = round(linspace(1, n_lags, 5));
% %tick_labels = arrayfun(@(x) sprintf('%.1fms', time_range(x)), tick_indices, 'UniformOutput', false);
Lag_Active = linspace(0, 400, n_lags);
Lag_obs = linspace(0, 400, n_lags);
% 
% %%% colori piu sfumati.....OBS--> ACT
sigma_gauss=1;
%HEATMAP ACT --> OBS
% Image Processing Toolbox
%M_smooth = imgaussfilt(corr_avg, sigma_gauss); 
imagesc(Lag_obs, Lag_Active, corr_avg);
% 
% %imagesc(Lag_obs, Lag_Active, corr_matrix);
set(gca, 'YDir', 'normal');  
colormap(redblue256());
%colormap(redbluecmap(256));   % Se non hai redbluecmap, leggi sotto per crearla

colorbar;
caxis([-1, 1]);  
y_lab=['Lag (ms) ',upper(config.active_subject)];
x_lab=['Lag (ms) ',upper(config.obs_subject)];
%y_lab=['Lag (ms) ',upper(config.obs_subject)];
%x_lab=['Lag (ms) ', upper(config.active_subject)];
% , ' - Obs'] o act
xlabel(x_lab);
ylabel(y_lab);
%title_ = sprintf('ACT-->OBS \nCondition %d: \nDirection %d',config.condition,  dir_);
title_ = sprintf('\nCondition %d: \nDirection %d',config.condition,  dir_);
%title_=['Condition ', num2str(config.condition), ': ' upper(config.active_subject), ... 
%    ' Act ',' - ' upper(config.obs_subject), ' Obs', 'Direction\n', num2str(dir_)]
title(title_);
% diagonal
hold on
xlims = xlim;          % [x_min  x_max]   (asse Obs)
ylims = ylim;          % [y_min  y_max]   (asse Active)
%plot(xlims, ylims, 'k-', 'LineWidth', 1.5);   % linea nera continua
% 
plot(xlims, ylims, '--', 'Color',[.5 .5 .5], 'LineWidth',1.2)
% 
hold off
%  output_dir = fullfile('HeatMaps', 'full_maps');
%  saveas(gcf, fullfile(output_dir, sprintf('act_to_obs_LagLag_condition%d_dir%d.png', config.condition, dir_)));
%  saveas(gcf, fullfile(output_dir, sprintf('AVG_LagLag_condition%d_dir%d.png', config.condition, dir_)));
 
%% Calcolo medie varie
corr_matrix_ = corr_avg;
n_lags = size(corr_matrix_, 1);

% Maschere per isolare le parti di interesse (sopra sotto e sulla diagonle
% principaòe
upper_triangle = triu(corr_matrix_, 1);       
lower_triangle = tril(corr_matrix_, -1);      
diagonal_values = diag(corr_matrix_);        

% Rimuovo valori che non vanno nel calcolo della medai
eps = 1e-6;
upper_triangle(abs(upper_triangle) < eps) = NaN;
lower_triangle(abs(lower_triangle) < eps) = NaN;
diagonal_values(abs(diagonal_values) < eps) = NaN;

%MEDIE GLOBALI
mean_total         = mean(corr_matrix_(:), 'omitnan')
mean_upper_global  = mean(lower_triangle(:), 'omitnan')
mean_lower_global  = mean(upper_triangle(:), 'omitnan')
mean_diag          = mean(diagonal_values, 'omitnan')

% MEDIE PER RIGA (asse Y / Active)
mean_upper_vec_rows = mean(lower_triangle, 2, 'omitnan'); 
mean_lower_vec_rows = mean(upper_triangle, 2, 'omitnan');

%MEDIE PER COLONNA (asse X / Obs)
mean_upper_vec_cols = mean(lower_triangle, 1, 'omitnan');  
mean_lower_vec_cols = mean(upper_triangle, 1, 'omitnan');
 
%   STORE
cond = config.condition;
results(cond).cond(dir_).mean_total            = mean_total;
results(cond).cond(dir_).mean_upper_global     = mean_upper_global;
results(cond).cond(dir_).mean_lower_global     = mean_lower_global;
results(cond).cond(dir_).mean_diag             = mean_diag;
results(cond).cond(dir_).mean_upper_vec_rows   = mean_upper_vec_rows;
results(cond).cond(dir_).mean_lower_vec_rows   = mean_lower_vec_rows;
results(cond).cond(dir_).mean_upper_vec_cols   = mean_upper_vec_cols;
results(cond).cond(dir_).mean_lower_vec_cols   = mean_lower_vec_cols;

%% Barplots tutti i risultati BARPLOTS 
results=load('results_29_09_25_all_trials.mat')
results=results.results
conds = 1:3;      
cc=numel(conds)
direzioni = 1:8;   
dd=numel(direzioni)
%results_18_08= results
mean_upper_all=zeros(cc,dd);
mean_lower_all=zeros(cc,dd);
mean_diag_all=zeros(cc,dd)
sd_upper_all=zeros(cc,dd)
sd_lower_all=zeros(cc,dd)
sd_diag_all=zeros(cc,dd)
% itero su condizioni e direzioni
for c = 1:cc
    %display (results(c))
    for d = 1:dd 
               elem = results(c).cond(d);
               display(results(conds(c)).cond)
               diffs(c,d) = elem.mean_K2S - elem.mean_S2K;
               diags(c,d) =elem.mean_diag;
               mean_upper_all(c,d ) = results(conds(c)).cond(d).mean_K2S;
               mean_lower_all(c, d) = results(conds(c)).cond(d).mean_S2K
               mean_diag_all(c, d)  = results(conds(c)).cond(d).mean_diag;
               sd_upper_all(c,d ) = results(conds(c)).cond(d).sd_K2S;
               sd_lower_all(c, d) = results(conds(c)).cond(d).sd_S2K;
               sd_diag_all(c, d) = results(conds(c)).cond(d).sd_diag

    end
 end

% Figure con 3 subplot (uno per upper/lower/diag)
figure;
subplot(3,1,1);
Y_u=mean_upper_all'
sd_u=sd_upper_all'
h1=bar(Y_u, 'grouped');
hold on 
for i = 1:numel(h1)                 
    x = h1(i).XEndPoints;           
    %errorbar(x, Y_u(:,i), sd_u(:,i), sd_u(:,i), ...
     %        'k', 'linestyle','none', 'LineWidth',1, 'CapSize',6);
end
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
for i = 1:numel(h2)
    x = h2(i).XEndPoints;
    %errorbar(x, Y_l(:,i), sd_l(:,i), sd_l(:,i), 'k','linestyle','none','LineWidth',1,'CapSize',6);
end
title('S → K');
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
for i = 1:numel(h3)
    x = h3(i).XEndPoints;
    %errorbar(x, Y_d(:,i), sd_d(:,i), sd_d(:,i), 'k','linestyle','none','LineWidth',1,'CapSize',6);
end
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
% conds = 1:3;        % condizioni
% direzioni = 1:8;    % direzioni
n_dir=numel(direzioni)
n_cond=numel(conds)
 % results_18_08= results
% Cotnainer per differenze e valori ion diagonale
diffs=zeros(n_dir,n_cond )
diags=zeros(n_dir,n_cond)
%itero su condizioni e direzioni
for d = 1:n_dir 
    for c = 1:n_cond 
        elem = results(c).cond(d);
        diffs(d,c) = elem.mean_K2S - elem.mean_S2K;
        diags(d,c) =elem.mean_diag;
    end
end

% 2) soglia
% dj base, se la differenza è inferiore a un certo valore, 
% non gli do peso (nel senso non è uno scostamento tale da poter dire
% in modo netto che un soggetto anticipa l'altro...if any)
epsThr = 0.05;                        
isNS   = abs(diffs) < epsThr;          
% 
% 3) etichette basate sul segno con soglia
% creo matrice di stringhe di dimensione condizioni per direzioni
signlab = strings(n_dir, n_cond);
signlab(isNS)                = "n.s.";
signlab(~isNS & diffs > 0)   = "K → S";
signlab(~isNS & diffs < 0)   = "S → K";
signlab=fliplr(signlab)
% % 
% % 4)layout plot orizzontale, un pannello per  direzione
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% ciclo si valori 1...n creando stringa dir i
% che vanno sul lato six per indicare le direZioni 
%condlabels = arrayfun(@(k) sprintf('cond %d', k), 1:n_cond, 'UniformOutput', false);
% Numero di condizioni
n_cond = 3;

% Etichette personalizzate per ogni condizione
labels = { ...
    'Cond 1: S act - K obs', ...
    'Cond 2: K act - S obs', ...
    'Cond 3: S act - K act' ...
};

% Ora creiamo le label in base a n_cond
condlabels = labels(1:n_cond);
condlabels = flip(condlabels);    

figure; clf;
% 2 colonne
ncols = 2;    
% righe necessarie
nrows = ceil(n_dir / ncols);       
tiledlayout(nrows, ncols, 'TileSpacing','compact', 'Padding','compact');
for c = 1:n_dir
    nexttile; hold on;
     % Inverto i dati per matchare l'ordine visuale (cond 1 in alto)
    data = diffs(c, :);
    data = flip(data);
    % barra
    h = barh(data, 'FaceAlpha', 0.95);
    title(sprintf('dir %d', c));

    % ASSE GRIGLIA TITOLO
    if mod(c, 2) == 1
        yticks(1:n_cond);
        yticklabels(condlabels);
    else
        yticks([]);   % rimuove le etichette Y per i subplot a destra
    end

    xlabel('Corr Diff');
    %yticks(1:n_cond); yticklabels(condlabels);
    xlim([-0.5 0.5]); grid on;

    % Colori: verde >0, rosso <0, grigio per n.s.
    cd = repmat([0.7 0.7 0.7], n_cond, 1);  % parte tutto grigio [n_cond x 3]
    idxPos = ~isNS(c,:) & diffs(c,:) > 0;  % barre significative positive
    idxNeg = ~isNS(c,:) & diffs(c,:) < 0;  % barre significative negative
    cd(idxPos,:) = repmat([0.2 0.6 0.2], sum(idxPos), 1);  % verde
    cd(idxNeg,:) = repmat([0.8 0.2 0.2], sum(idxNeg), 1);  % rosso
    h.FaceColor = 'flat';      
   % h.CData     = cd;         
    h.CData     = flip(cd,1);

    % Etichette a DESTRA 
    for d = n_cond:-1:1
        % posizione a destra del grafico, indipendente dal valore
        % aggiunge uno 0.qcs limiti
        xpos = max(xlim) + 0.01;  
        %scrive lsetichetta sulla riga d, allineata in mezzo verticalmente,
        % e “appoggiata” a sinistra orizzontalmente (verso l%interno del
        % grafico).
        text(xpos, d, signlab(c,d), ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left', ...
        'FontWeight', 'bold');
    end

    hold off;
end
% 
% 
%%%% PLOTTIAMO i valori di correlazione sulla diagonale
figure; clf;
% 2 colonne
ncols = 2;    
% righe necessarie
nrows = ceil(n_dir / ncols);       
tiledlayout(nrows, ncols, 'TileSpacing','compact', 'Padding','compact');
for c = 1:n_dir
    nexttile; hold on;
     % Inverto i dati per matchare l'ordine visuale (cond 1 in alto)
    data = diags(c, :);
    data = flip(data);
    % barra
    h = barh(data, 'FaceAlpha', 0.95);
    title(sprintf('dir %d', c));

    % ASSE GRIGLIA TITOLO
    if mod(c, 2) == 1
        yticks(1:n_cond);
        yticklabels(condlabels);
    else
        yticks([]);   % rimuove le etichette Y per i subplot a destra
    end

    xlabel('Corr Diff');
    %yticks(1:n_cond); yticklabels(condlabels);
    xlim([-0.5 0.5]); grid on;
end





