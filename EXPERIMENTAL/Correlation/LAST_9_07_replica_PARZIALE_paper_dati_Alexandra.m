
%% Ref libreria python https://zenodo.org/records/12359299
% parametrizzare blocchi
% label temporale
% distinzione 

% quando scrivo "per loro", intendo per Hasson e compagnia cantante

% Replica moooolto semplificata del paper
% mi sono focalizzato sui punti chiave; 
% 1) usiamo i dati s e k delle scimmie nelle reciproche condizioni act e obs;
%    tutta la parte di allineamento dei dati (dialoghi e dati neurali nel
%    loro caso, movimenti del braccio e dati neurali nel nostro, sono
%    saltati.
% 2) I dati sono ricampionati e ridotti a un quinto della lunghezza
%     orginaria per avere dei dati sensati/decenti
%     (formula K=[(N-W)/S]+1...per ogni trial
%     K Numero di punti rimanenti (118 obs (una per 5 ms) per un totale di 599 ms)
%     N lunghezza dati originari (599 obs (una per ms) obs per trial), 
%     W(indow) ampiezza finestra (10), 
%     S(hift) passo di avanzamento (5) 
%     

% 3)  Genero gli embeddidng con CEBRA behaviour (diciamo supervised)
%     usando più o meno sempre la stessa struttura: 
%     batch_size: 1024
%     max_iterations: 8000-12000
%     hybrid: False
%     verbose: True
%     time_offsets: 10
%     output_dimension: 3
%     learning_rate: 0.0001
%     num_hidden_units: 32-64
%     temperature: 1

% Gli embedding sono il corrispettivo di quello che loro generano con le
% parole usando gpt2. Dati proiettati su uno spazio diverso.
% Invero nel loro caso è diverso dacchè lo spazio generato è esterno a
% soggetti
% Il fatto che ci siano parole (movimenti) ripetuti non rileva in quanto si
% tratta (in teoria) di manifold contestuali

% 2) ho allineato i dati campionando come lo fanno loro
% 3) regressione ridge per lag e per elettrodo con cross validation esterna 
% e nested per ricerca lambda ottimale
%    e calcolo correlazioni tra valori predetti e osservati -(correlazione
%    tra segmenti dello stimolo). Poco senso con i dati che abbiamo. 
% 4) Brain to brain coupling 
%  praticamente testano il modello allenato su un soggetto, sull'altro
%  soggetto...

% Le parti di riferimento sul paper sono nei metodi del paper:
% "A shared model-based linguistic space for transmitting our thoughts from 
%          brain to brain in natural conversations"
% Encoding analysis
% Model-based brain-to-brain coupling

%% creo struttura per accesso dinamico ai dati 
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
% QUESTE SONO LE INFO INIZIALI Mutevoli  
% numero trial e conseguente  Numero di labels 
% 64 trials (totali, 8 per direzione)  da 599 ms (e una osservazione per ms)
% trials per direzioni
%n_label=8 ;
n_trials=64 ; 
% lunghezza trial originaria 
trial_length=599;
% ricampionati con finestra  variabile (parametri ricampionamento)
% finestra di ricampionamento: 10 ms. Passo: 5 
% valori usati per resampling (se metti tutto a 1, ricadi nella lunghezza 
% originaria dei dati)
shift_=5 ;
window_=10 ;
% lunghezza trial resamplati
trial_len_res=floor((trial_length-window_)/shift_)+1;

% durata sessione in secondi (64 trials da 599 ms/bin...ridotti a 118 bin)
% circa 38" totali (se considero tutti i trial nella loro interezza)
% per loro erano 60" di dialogo
dd_total= n_trials*trial_length/1000; 
% frequenza campionamento (numero di punti al secondo...dopo aver
% ricampionato abbiamo un quinto dei dati circa)
%  % per loro in hertz sono 512 (punti al secondo)
fs=ceil(1000/(trial_length/trial_len_res));

% rimuovere  tutti i dati prima dello zero per ogni trial
new_indices = [];
zeros_=200;
% vado per eccesso
zero_index_resampled = ceil((zeros_ / 599) * trial_len_res);

for t = 0:(n_trials - 1)
        start_idx = t * trial_len_res + 1;
        trial_indices = start_idx : (start_idx + trial_len_res - 1)
% solo la parte post-zero (dalla 40ª bin in poi)
        new_indices = [new_indices, trial_indices(zero_index_resampled:end)];
end

%% Carichiamo dati relativi a una specifica condizione
%(una scimmia agisce. l'altra osserva)
% condizione 1: S agisce K osserva
% Condizione 2: K agisce S osserva
% Condizione 3: S agisce K agisce
% condizioni > 3 varianti eventuali: 
% cond 4: sono tutte e 3 le condzioni precedenti messe insieme (da inserire)
%%%%%%%%%%%%%%%%%%%%%%%%% LOAD DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

config.condition =  1;
config.active_subject = 's';
config.obs_subject = 'k';

[active_neural, active_label, active_embed, ...
 obs_neural, obs_label, obs_embed] = get_data_from_config(config, DATA, new_indices);

%%%%%%%%%%%%%%%%%%%%%% filtra per direzione %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  filtro pr direzione (1-8), cambiano una serie di cose
dir_=2;
n_label=1;
%%% 8 trials per direzione  (per trial indico il singolo gesto)
% 24 quando mettiamo inismee i dati (da fare)
n_trial=8 
%%  lunghezza trial in tempo
trial_length_=trial_length-zeros_+1
% filtrare poi i dati partedno dalle label e rispettive pos
label_indices=find(active_label==dir_)
% Lunghezza in bin resamplatei del singolo trial partendo da zero
trial_len_res_=trial_len_res-zero_index_resampled+1
% tempo per sessione (sessione = 8 trials per direzione
%  un singolo trial
tt_single = (0:trial_len_res_-1) * shift_ + window_;  
tt = zeros(1, n_trial * trial_len_res_);
for i = 1:n_trial
    start_idx = (i-1) * trial_len_res_ + 1;
    stop_idx  = i * trial_len_res_;

    trial_start = (i-1) * trial_length_;  % in ms

    tt(start_idx:stop_idx) = tt_single + trial_start;
end     
% tempo in secondi della sessione (sempre 8 trials per driezione)
% 3.2 secondi
tt = tt / 1000;  % 


% Dati neruali di partenza
% % act: la scimmia  che muove il braccio
% % obs:e' la scimmia che osserva
active_neural=active_neural(label_indices,:);
obs_neural=obs_neural(label_indices,:);
active_label=active_label(label_indices,:);
obs_label=obs_label(label_indices,:);
% manifold genrate con CEBRA (o altro)
active_embed = active_embed(label_indices,:);
obs_embed = obs_embed(label_indices,:);


% embeddings MANIFOLD
% - per loro sono generati con GPT2 (che per ogni parola tira fuori un 
% embeddding di 1600 componenti)...
% per noi sono quelli di cebra qui metto meno componenti

% numero elettrodi....da capire come usare e adattare questa info poichè 
% nel mio  caso non sono lo stesso numero (per ora li ho scelti con shap)
%n_channels=;

% onset dei movimenti (nel loro caso si considera un intorno di durata arbitraria
% tipo -+4 sec per ogni parola lungo la durata complessiva del dialogo) 
% nel nostro setting abbiamo un gesto che inizia (in ogni trial)
% 100 ms dal principio del trial. Inoltre, invece che considerare un
% intorno del movimento consideriamo solo momenti successivi (ref. Fra)
% considerando 64 trial di 118 momenti (5ms) so che il movimento inizia
% alla fine della 19ma obs per trial (91-100 ms).

%intorno_sec=4;
%word_onsets = sort(randi([intorno_sec/2*fs length(tt)-intorno_sec/2*fs], 1, wpm));

%%% devo trovare la finestra dove cade l'onset del movimento (per trial)
% Con i dati resamplati  è la prima per trial. 
onset_bin=find(tt>=0,1)
%%% 
mov_onsets=tt(onset_bin:trial_len_res_:end)


%% ALLINEAMENTO DATI NEURALI - DATI movimento Embedding %%%%%%%%%%%%%%%%
% di fatto si analizza l'attività neurale in intervalli di tempo specifici
% successivi al principio del movimento (Hasson&Co. lo fanno intorno
% agli onsets delle parole - 4 secondi prima e 4 dopo ogni parola).

% come loro, che creano finestre di 250ms con overlap di 62.5 ms, creiamo
% finestre con overlap; in particolare,  considero finestre di block_size ms 
%  con overlap di block_size-shift ms 


% Analogamente procediamo per manifold generate con CEBRA.
% questo è altro elemnto distintivo dal paper in quanto loro hanno un
% embedding comune generato con GPT (oltretutto su uno spazio molto più
% grande) mentre noi abbiamo 2 manifold (1 per scimmia) su 3 dimensioni

% inoltre, il loro embedding è statico...il nostro si muove nel tempo; per
% loro ogni parola ha un suo embedding multidimensionale nello spazio ma
% unidimensionale nel tempo; i nostri embedding si muovono nel tempo e
% nello spazio (per avere una situazione analoga alla loro dovrei avere per
% ogni movimento una solo osservazione nello spazio 3d). Da un lato questo
% torna comodo in quanto alla fine mi ritrovo con strutture analoghe per
% dati neurali ed embedding (per ricondurci al loro caso, anche sulle
% manifold generate con cebra miedamo sul tempo...più avanti)

% sono bin da 5 ms 
% da cui se block size è bb ms e shift block ss ms, sono blocchi con 
% overlap di bb-ss ms 
% 1-bbms, 1+ss-bb+ss_ms e così.....

% la scelta delle dimmensioni dei blocchi (lags)è del tutto arbitraria...
% cerchiamo di includere tutti i dati. 
% Se non includo lo zero (quidni 5ms dopo stante il resampling) (18-5 ad esempio)
% Qui parto includendo lo zero
block_size=19
shift_block=6

for i = 1:n_trial
    onset_pos = find(abs(tt - mov_onsets(i)) < 1e-6, 1);
    
    % se non includiamo lo zero +1 qui sotto
    idx_min =onset_pos;
    % 
    idx_max = onset_pos+trial_len_res_-onset_bin;
    start_idx = idx_min;
    blocchi_idx = [];
    
    while start_idx + block_size - 1 <= idx_max
        end_idx = start_idx + block_size - 1;
        blocchi_idx = [blocchi_idx; start_idx, end_idx]
        
        % Calcola il prossimo indice di partenza con shift
        start_idx = start_idx + shift_block;
    end
      
    dati_obs{i} = cell(size(blocchi_idx, 1), 1);
    dati_active{i} = cell(size(blocchi_idx, 1), 1);
    manif_active_match{i} = cell(size(blocchi_idx, 1), 1);
    manif_obs_match{i} = cell(size(blocchi_idx, 1), 1);
    for j = 1:size(blocchi_idx, 1)
        start_blocco = blocchi_idx(j, 1);
        end_blocco = blocchi_idx(j, 2);
        dati_active{i}{j} = active_neural(start_blocco:end_blocco, :);
        dati_obs{i}{j} = obs_neural(start_blocco:end_blocco, :);
        manif_active_match{i}{j}=active_embed(start_blocco:end_blocco, :);
        manif_obs_match{i}{j}=obs_embed(start_blocco:end_blocco,:);
    end
end

% quidni per ogni scimmia mi ritrovo con strutture (dati obs e dati active
% di dimensione 1 (soggetto) per 8 (trials per direzione) e dentro ogni 
% cella ho i dati neurali con riferimento il principio del movimento:
% quindi 8 movimenti...ogni movimento ha un tot blocchi di dimensione block
% size.  Ogni blocco rappresenta quindi osservazioni neurali "laggate" nel
% tempo

% per quanto riguarda le manifold anche qui, 8 celle (per direzione)
% che contengono i dati cebra momenti contemporanei ai dati neurali estratti

%% RIDGE REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nel paper originale, l'embedding (X) è statico per ogni parola: 
% dimensione [#parole × spazio_embedding GPT2 ad esempio 500]
% L'attività neurale (Y) è laggata: per ogni lag e per ogni neurone 
% si fa la media sul tempo (durata) del blocco (lag) neurale.
% Il risultato è un vettore Y di dimensione [parole × 1] per ogni combinazione
% (neurone × lag). Ad esempio ipotizzando 150 parole, 60 blocchi temporali,
% o lags, intorno all'onset di ogni parola e 64 neuroni/canali, avrò 60*64 y
% di  dimensione 150*1. La X sarà 150*dimensione embedding
% Si esegue quindi una regressione per ciascun neurone 
% e per ciascun lag,  usando sempre la stessa matrice X.

% Nel nostro caso,  L'attività neurale (Y) ha la stessa struttura: 
% una cella per movimento (trial),; ogni cella contiene blocchi/lag di 
% dimensione  [t_lag × #neuroni] N.B.! T_lag=block_size!
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
% (saranno tanti vettori di dimensione (numero movimenti*1) 
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

for b = 1:n_blocchi
    for e = 1:n_channels_active
        Y_active{b,e} = [];
        %Y_l{b,e} = [];
        for p = 1:n_trial         
            Y_active{b,e} = [Y_active{b,e}; mean(dati_active{p}{b}(:, e))];
            %Y_l{b,e} = [Y_l{b,e}; mean(blocchi_dati{p}{b}(:, e))];
        end
    end
end

% X_active
for b = 1:n_blocchi
        X_active{b} = [];
        %Y_l{b,e} = [];
        for p = 1:n_trial
            
            X_active{b,1} = [X_active{b,1}; mean(manif_active_match{p}{b})];
            %Y_l{b,e} = [Y_l{b,e}; mean(blocchi_dati{p}{b}(:, e))];
        end
 end

%%%%%%%%%%%%%%% SCIMMIA passiva
Y_obs= cell(n_blocchi, n_channels_obs);
X_obs= cell(n_blocchi, 1);

% Y_s Riempimento di Y_S con i dati
for b = 1:n_blocchi
    for e = 1:n_channels_obs
        Y_obs{b,e} = [];
        %Y_l{b,e} = [];
        for p = 1:n_trial
            Y_obs{b,e} = [Y_obs{b,e}; mean(dati_obs{p}{b}(:, e))];
        end
    end
end

% X_obs
for b = 1:n_blocchi
    
        X_obs{b} = [];
        %Y_l{b,e} = [];
        for p = 1:n_trial
            
            X_obs{b,1} = [X_obs{b,1}; mean(manif_obs_match{p}{b})];
        end
    end

% Alla fine per ogni scimmia mi ritrovo con due strutture Y e X. La Y
% contiene tanti vettori di dimensione #trials*1 (se mantenessi il tempo 
% ogni vettore sarebbe espanso in formato long e avrebbe una dimensione 
% t_lag*#movimenti*1)

% % In particolare 
% facciamo una regressione per lag e per neurone come nel paper originale
% la X contiene le manifold allineate nel tempo (anche qui ho fatto media 
% sulla lungheza del blocco) per blocco temporale. 

%% Correlation analysis
% Questo di fatto è il punto 1 nel paragrafo Model-based brain-to-brain
% coupling  MEtodi (pag e4)
% "we develop a framework for assessing five types of generalization simultaneously:
% testing encoding model generalization
% (1) across segments of the stimulus (using 10-fold cross-validation)"
% Ora, in teoria, il modello viene testato utilizzando una 10 fold
% cross-validation 
% In ogni iterazione, un sottoinsieme di blocchi di parole viene escluso 
% dal training e utilizzato come test set 
% Dopo l'addestramento Le risposte neurali predette (`y_hat`) per il test
% set vengono calcolate moltiplicando la matrice delle embedding del test 
% set con i coefficienti stimati. Si calcola la correlazione tra i valori
% predetti y_hat e le osservate.
% Il processo viene ripetuto per tutti i fold della cross-validation.
% Alla fine, si ottiene una correlazione media tra predetto e osservato per
% ogni fold,che rappresenta la performance del modello nell'encoding delle
% risposte neurali

% Io faccio una cosa analoga con una cross validation nested per la scelta 
% lambda ottimale GCV per ridurre la complessità computazionale e
% soprattutto per usare tutti i dati ad ogni lag.
% la GCV viene fatta su dati training per ogni fold

%%
% RIDGE REGRESSION per stimare i parametri utili a B2B coupling 
% NO Cross validation (solo GCV per la scelta del lambda ottimo)
% ref per GCV: "The Elements of Statistical Learning " (Hatie et, al 2008) 
%          Generalized Cross-Validation as a Method for Choosing a Good
%          Ridge Parameter, (Golub et al., 1979 Technometrics)

% y length for each electrode-lag combo
t_lag=length((dati_obs{1}{1}(:, 1)));
n_movements=n_trial;

% Inizializzazione delle strutture per i coefficienti, predizioni e
% correlazioni  
% Scimmia active
coefs_active = cell(n_blocchi, n_channels_active);
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
        y_hat = X_ * beta;
        %%% da creare contenitore prima
        corrs_active(b, e) = compute_metric(y_, y_hat, y_);
           
    end 
end


% come loro prendo una soglia (che loro caclolano empiricamente con permutazioni
% e tengo solo le correlazioni oltre quella soglia--
r_threshold = 0.05;  
corr_active_ = abs(corrs_active) > r_threshold;
corr_active_ = corrs_active .* corr_active_ ;

%%%%%%%%%%%%  Scimmia obs
coefs_obs =cell(n_blocchi, n_channels_obs);
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
        y_hat = X_ * beta;
        %%% da creare contenitore prima
        corrs_obs(b, e) = compute_metric(y_, y_hat, y_);           
    end 
end

% come loro prendo una soglia (che loro caclolano empiricamente con permutazioni
% e tengo solo le correlazioni oltre quella soglia--
r_threshold = 0.05;  
corr_obs_=  abs(corrs_obs) > r_threshold;
corr_obs_ = corrs_obs .* corr_obs_ ;

%%
%%%%%%%%%%%%%%%%%%%%%%% 4) Brain  to Brain coupling %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% punto 2 nel paragrafo Model-based brain-to-brain coupling MEtodi (pag e4)
%%%%%%%%%%%%%%%%%%  GENERALIZATION ACROSS SUBJECTS ######################
% "we develop a framework for assessing five types of generalization simultaneously:
% testing encoding model generalization
% (2) across subjects (within speaker–listener dyads)"
% Generalizzazione inter-soggetto (cross-subject):
% Uso i coefficienti stimati su un soggetto per predire l'attività neurale
% dell'altro, sfruttando i rispettivi embedding temporali (X). 
% A differenza del within-subject, qui le X (manifold) sono diverse
% tra attivo e passivo, ma sincronizzate nel tempo.
% Questo permette di testare se il pattern neurale è condiviso tra i soggetti
% (ad es. act ricostruisce i pattern neurali di obs e viceversa).

% (da rifare passando per SHAP su entrambi i soggetti con ordinamento
% quindi stabilito da SHAP)
% Per questa parte mi occorrono gli stessi neuroni, o meglio, un numero di
% neuroni (canali) uguale per scimmia...dal momento che K ha 202 neuroni e
% S ne ha 184, seleziono 184 nuroni per K. Lo faccio sulla base degli shap
% values (in python) con il quale ho creato un "ranking" di importana di
% sulla base del quale ne tolgo 18

% Neuroni da rimuovere (la scimmia  s ha piu neuroni di k)
to_remove=[97,156,109,165,25,58,153,4,190,10,33,65,3,135,103,48,118,186]
if n_channels_obs>n_channels_active
    X_active_r=X_active
    Y_active_r=Y_active
    X_obs_r=X_obs
    Y_obs_r=Y_obs
    % Protezione da indici fuori range
    to_remove = to_remove(to_remove <= size(Y_obs_r, 2));
    Y_obs_r(:,to_remove)=[]
    coefs_active_r=coefs_active
    coefs_obs_r=coefs_obs
    coefs_obs_r(:, to_remove,:)=[]
else
    X_obs_r=X_obs
    Y_obs_r=Y_obs
    X_active_r=X_active
    Y_active_r=Y_active
    % Protezione da indici fuori range
    to_remove = to_remove(to_remove <= size(Y_active_r, 2))
    Y_active_r(:,to_remove)=[]  
    coefs_obs_r=coefs_obs
    coefs_active_r=coefs_active
    coefs_active_r(:, to_remove,:)=[]

end

n_channels=length(Y_active_r)

%% Si parte da strutture blocchi (lags) * neuroni contenenti ognuna 
% i coefficienti stimati within subject. Si usano questi coefficienti 
% reciprocamente per proiettare l'attività neurale di un soggetto nello
% spazio degli embedding dell'altro

corrs_inter_active_obs = zeros(n_blocchi, n_channels);
corrs_inter_obs_active = zeros(n_blocchi, n_channels);

for b = 1:n_blocchi
    X_active_ = X_active_r{b};
    X_obs_    = X_obs_r{b};

    for e = 1:n_channels
        y_active_ = Y_active_r{b, e};
        y_obs_    = Y_obs_r{b, e};

        % Uso i beta mediati invece che quelli fold‐specifici
        % usato per predire active
        beta_obs   = coefs_obs_r{b, e};  
        % usato per predire obs
        beta_active = coefs_active_r{b, e};  

        if isempty(coefs_obs_r) || isempty(coefs_active_r)
            continue
        end

        % Predizione (tutti i dati, no test/train perché già validato)
        y_hat_active = X_active_ * beta_obs;
        y_hat_obs    = X_obs_ * beta_active;

        % Correlazioni tra predetto e osservato (inter-subject)
        % quello che otteniamo sono delle matrici lag*chns in cui ogni riga
        % i rappresenta la correlazione media per tutti i canali al lag i
        % ccioè quanto il soggetto obs (il modello addestrato su obs invero
        % è in grado di predire il soggeto active con shift temporale i 
        % (e viceversa 
        corrs_inter_active_obs(b, e) = compute_metric(y_active_, y_hat_active);
        corrs_inter_obs_active(b, e) = compute_metric(y_obs_, y_hat_obs);
    end
end

%%
%%%%%% COSTRUZIONE HEATMAPS lag lag fig 3a %%%%%%%%%%
n_lags=n_blocchi, n_channels


% N.B.!!  media sul canale? o sulle fold (in articolo origianel lavorano
%  su fold che noi non abbiamo data la poca disponibilità di dati
% quidni noi ci ritroviamo con matrici di correlazione lags*canali

% Soglia minima per eliminare zeri o valori troppo piccoli
epsilon = 1e-4;

% Pulizia di corrs_inter_active_obs 
mask_small_active_obs = abs(corrs_inter_active_obs) < epsilon;
corrs_inter_active_obs(mask_small_active_obs) = ...
    sign(corrs_inter_active_obs(mask_small_active_obs)) * epsilon;

% Pulizia di corrs_inter_obs_active 
mask_small_obs_active = abs(corrs_inter_obs_active) < epsilon;
corrs_inter_obs_active(mask_small_obs_active) = ...
    sign(corrs_inter_obs_active(mask_small_obs_active)) * epsilon;

% Inizializzo la matrice lag-lag
corr_matrix = zeros(n_lags, n_lags);

% calcolo correlazione con fold
% Calcolo  l correlazione per ogni combinazione di lag K-S
%  l'indice i rappresenta il lag usato per il flusso active→ obs
%  da cui, l'indice j rappresenta il lag usato per il flusso obs→active
for i = 1:n_lags
    for j = 1:n_lags
        % Vettore di valori per il lag i (Active to obs)
        % (1 × folds)
        Active_vector = corrs_inter_active_obs(i, :); 
        
        % Vettore di  valori per il lag j (obs to active)
        % (1 × folds)
        Obs_vector = corrs_inter_obs_active(j, :); 
        
        % Calcolo la correlazione tra i due vettori
        %Correlazione fra il vettore dei valori (calcolato per 
        % Active→Obs al lag i) e il vettore per Obs→Active al lag j.
        %   lag‐lag matrix (n_lag x n_lags).
        corr_matrix(i, j) = corr(Active_vector', Obs_vector');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT PLOT FULL %%%%%%

% time range (vado a spanne...da perfezionare...sono circa 400 ms dopo il movimento)
% Array dei tempi associati ai lag
%time_range = linspace(0, 400, n_lags);
% Seleziona 5 tick equidistanti
%tick_indices = round(linspace(1, n_lags, 5));
% %tick_labels = arrayfun(@(x) sprintf('%.1fms', time_range(x)), tick_indices, 'UniformOutput', false);
Lag_Active = linspace(0, 400, n_lags);
Lag_obs = linspace(0, 400, n_lags);
% 
% %%% colori piu sfumati
sigma=1;
M_smooth = imgaussfilt(corr_matrix, sigma);  % Image Processing Toolbox
imagesc(Lag_obs, Lag_Active, M_smooth);
% 
% %imagesc(Lag_obs, Lag_Active, corr_matrix);
set(gca, 'YDir', 'normal');  
%colormap('jet');
colormap(redbluecmap);   % Se non hai redbluecmap, leggi sotto per crearla

colorbar;
caxis([-0.25, 0.25]);  
y_lab=['Lag (ms) ',upper(config.active_subject),' - Act'];
x_lab=['Lag (ms) ',upper(config.obs_subject),' - Act'];
xlabel(x_lab);
ylabel(y_lab);
% %%% sprint f per interpretare a capo e stringa numeri intere
title_ = sprintf('Condition %d: \nDirection %d',config.condition,  dir_);
% %title_=['Condition ', num2str(config.condition), ': ' upper(config.active_subject), ... 
%  %    ' Act ',' - ' upper(config.obs_subject), ' Obs', 'Direction\n', num2str(dir_)]
title(title_);
% %%% diagonal
hold on
xlims = xlim;          % [x_min  x_max]   (asse Obs)
ylims = ylim;          % [y_min  y_max]   (asse Active)
%plot(xlims, ylims, 'k-', 'LineWidth', 1.5);   % linea nera continua
% 
plot(xlims, ylims, '--', 'Color',[.5 .5 .5], 'LineWidth',1.2)
% 
hold off

%%%% Salvataggio
output_dir = fullfile('HeatMaps', 'full_maps');
saveas(gcf, fullfile(output_dir, sprintf('LagLag_condition%d_dir%d.png', config.condition, dir_)));
% 

%% medie per i bar plot e time plot
%calcolare le correlazioni medie sopra e sotto la diagonale
corr_matrix_=corr_matrix


n_lags = size(corr_matrix_, 1);

% MASCHERE PER SOPRA/SOTTO ANTIDIAGONALE
% Fliplr per lavorare rispetto alla diagonale secondaria (da basso-sx ad alto-dx)
flipped = fliplr(corr_matrix_);

% Parte sopra antidiagonale escludendo diagonale
upper_flipped = triu(flipped, 1);         
upper_triangle = fliplr(upper_flipped); 

% Parte sotto antidiagonale
lower_flipped = tril(flipped, -1);
lower_triangle = fliplr(lower_flipped);

% Maschero per la diagonale secondaria (per escluderla anche a mano)
diag_idx = sub2ind([n_lags, n_lags], 1:n_lags, n_lags:-1:1);
diag_mask = false(n_lags);
diag_mask(diag_idx) = true;

% Imposto a NaN i valori sulla diagonale secondaria (antidiagonale)
upper_triangle(diag_mask) = NaN;
lower_triangle(diag_mask) = NaN;

%  MEDIE GLOBALI
mean_total             = mean(corr_matrix_(:), 'omitnan');
mean_upper_global      = mean(upper_triangle(:), 'omitnan');
mean_lower_global      = mean(lower_triangle(:), 'omitnan');

% MEDIE PER RIGA (centrate su ACTIVE lag) 
 % media riga per riga
mean_upper_vec_rows = mean(upper_triangle, 2, 'omitnan'); 
mean_lower_vec_rows = mean(lower_triangle, 2, 'omitnan');

% MEDIE PER COLONNA (centrate su OBS lag) 
% media colonna per colonna
mean_upper_vec_cols = mean(upper_triangle, 1, 'omitnan');  
mean_lower_vec_cols = mean(lower_triangle, 1, 'omitnan');

% Salvataggio in una struttura RESULTS
cond = config.condition;
results(cond).dir(dir_).mean_total            = mean_total;
results(cond).dir(dir_).mean_upper_global     = mean_upper_global;
results(cond).dir(dir_).mean_lower_global     = mean_lower_global;
results(cond).dir(dir_).mean_upper_vec_rows   = mean_upper_vec_rows;
results(cond).dir(dir_).mean_lower_vec_rows   = mean_lower_vec_rows;
results(cond).dir(dir_).mean_upper_vec_cols   = mean_upper_vec_cols;
results(cond).dir(dir_).mean_lower_vec_cols   = mean_lower_vec_cols;

%% per il plot troncato
% trunc_start= 4;
% trunc_end=4;
% n_lags_trunc=n_lags-trunc_start-trunc_end
% corr_matrix_trunc=corr_matrix(trunc_start+1:n_lags-trunc_end,trunc_start+1:n_lags-trunc_end)
% %% da migliorare e magari mettere insieme alcuni lag
% % time range (vado a spanne...da perfezionare...sono circa 400 ms dopo il movimento)
% % Array dei tempi associati ai lag
% 
% 
% tick_labels = {'100', '150', '200', '250', '300'};
% tick_pos = linspace(1, n_lags_trunc, 5);  % es. [1 7 13 19 25] se n_lags_trunc=25
% Lag_Active_trunc = linspace(100, 300, n_lags);
% Lag_obs_trunc = linspace(100, 300, n_lags);
% % Posizioni dei tick all’interno di n_lags_trunc (equidistanti tra 1 e n_lags_trunc)
% % Seleziona 5 tick equidistanti
% 
% %%% colori piu sfumati
% sigma=1
% M_smooth = imgaussfilt(corr_matrix_trunc, sigma);  % Image Processing Toolbox
% imagesc(M_smooth);
% 
% %imagesc(Lag_obs, Lag_Active, corr_matrix);
% set(gca, 'YDir', 'normal');  
% colormap('jet');
% colorbar;
% y_lab=['Lag (ms) ',upper(config.active_subject)];
% x_lab=['Lag (ms) ',upper(config.obs_subject)];
% xlabel(x_lab)
% ylabel(y_lab);
% 
% %  
% set(gca, 'XTick', tick_pos);  
% set(gca, 'XTickLabel', tick_labels);
% set(gca, 'YTick', tick_pos);
% set(gca, 'YTickLabel', tick_labels);
% 
% %xlabel(['Lag  (ms) ', upper(config.obs_subject)]);
% %ylabel(['Lag  (ms) ', upper(config.active_subject)]);
% %%% sprint f per interpretare a capo e stringa numeri intere
% title_ = sprintf('Condition %d: \nDirection %d', ...
%     config.condition,  dir_);
% %title_=['Condition ', num2str(config.condition), ': ' upper(config.active_subject), ... 
%  %    ' Act ',' - ' upper(config.obs_subject), ' Obs', 'Direction\n', num2str(dir_)]
% title(title_);
% %%% diagonale
% hold on
% xlims = xlim;          % [x_min  x_max]   (asse Obs)
% ylims = ylim;          % [y_min  y_max]   (asse Active)
% %plot(xlims, ylims, 'k-', 'LineWidth', 1.5);   % linea nera continua
% 
%  plot(xlims, ylims, '--', 'Color',[.5 .5 .5], 'LineWidth',1.2)
% 
% hold off





% 
% 
% n_lags_=n_lags
% corr_matrix_=corr_matrix
% 
% % --- Medie globali e per lag---
% % Sopra la diagonale
% % Sopra la diagonale
% upper_triangle = fliplr(triu(fliplr(corr_matrix_),-0));
% idx = sub2ind([n_lags_ n_lags], 1:n_lags, n_lags:-1:1);
% upper_triangle(upper_triangle==0)=NaN
% upper_triangle(idx) = NaN;
% 
% mean_upper_lag_OBS= mean(upper_triangle, 'omitnan');
% mean_upper_global_OBS=mean(upper_triangle(upper_triangle ~= 0), 'omitnan');
% mean_upper_lag_ACT= mean(upper_triangle,2, 'omitnan');
% mean_upper_global_ACT=mean(upper_triangle(upper_triangle ~= 0),'omitnan');
% 
% %mean_g= mean(mean_upper_lag, 'omitnan');
% 
% % Sotto la diagonale
% lower_triangle = fliplr(tril(fliplr(corr_matrix_), 0));
% idx = sub2ind([n_lags_ n_lags], 1:n_lags, n_lags:-1:1);
% lower_triangle(lower_triangle==0)=NaN
% lower_triangle(idx) = NaN;
% 
% mean_lower_lag_OBS= mean(lower_triangle, 'omitnan');
% mean_lower_global_OBS=mean(lower_triangle(lower_triangle ~= 0), 'omitnan');
% mean_lower_lag_ACT= mean(lower_triangle,2, 'omitnan');
% mean_lower_global_ACT=mean(lower_triangle(lower_triangle ~= 0),'omitnan');
% 
% % Tutta la matrice
% mean_total = mean(corr_matrix_(:), 'omitnan');
% % 
% % % ==== MEDIE PER RIGA (Active-lag centrico) ====
% % mean_upper_vec_rows = nan(n_lags_, 1);
% % mean_lower_vec_rows = nan(n_lags_, 1);
% % 
% % 
% % for i = 1:n_lags
% %     if i < n_lags
% %         mean_upper_vec_rows(i) = mean(corr_matrix_(i, i+1:end), 'omitnan');
% %     end
% %     if i > 1
% %         mean_lower_vec_rows(i) = mean(corr_matrix_(i, 1:i-1), 'omitnan');
% %     end
% % end
% % 
% % % ==== MEDIE PER COLONNA (Obs-lag centrico) ====
% % mean_upper_vec_cols = nan(n_lags_, 1);
% % mean_lower_vec_cols = nan(n_lags_, 1);
% % 
% % for j = 1:n_lags
% %     if j > 1
% %         mean_upper_vec_cols(j) = mean(corr_matrix_(1:j-1, j), 'omitnan');
% %     end
% %     if j < n_lags
% %         mean_lower_vec_cols(j) = mean(corr_matrix_(j+1:end, j), 'omitnan');
% %     end
% % end
% 
% cond=config.condition
% results(cond).dir(dir_).mean_total            = mean_total;
% results(cond).dir(dir_).mean_upper            = mean_upper;
% results(cond).dir(dir_).mean_lower            = mean_lower;
% 
% results(cond).dir(dir_).mean_upper_vec_rows   = mean_upper_vec_rows;
% results(cond).dir(dir_).mean_lower_vec_rows   = mean_lower_vec_rows;
% 
% results(cond).dir(dir_).mean_upper_vec_cols   = mean_upper_vec_cols;
% results(cond).dir(dir_).mean_lower_vec_cols   = mean_lower_vec_cols;
% 
% 
% % Allora effettivamente la diagonale è da dx a six (controintuitivo ma è
% % così)...ora andando nel dettaglio della matrice:
% % Correlazioni negative (in alto a destra), %Indicano che quando il segnale
% % active è già terminato o in calo (cioè a % lag elevato), il segnale 
% % obs continua invece ad essere alto. 
% % Dunque, il segnale obs "sopravvive" più a lungo rispetto ad active.
% % Questo si traduce matematicamente in una correlazione negativa, 
% % perché valori bassi (finali) di active sono associati sistematicamente a
% % valori ancora elevati di obs.
% %Correlazioni positive (in basso a sinistra):
% %Indicano invece che quando il segnale active è alto (o inizia a crescere), 
% % anche il segnale obs aumenta subito dopo. Questo significa che active 
% % anticipa obs.
% %La correlazione positiva è alta perché aumenti in active precedono
% % chiaramente gli aumenti in obs (o sono sincronizzati).
% 
% % Questa lettura è molto efficace e comunica chiaramente il rapporto 
% % temporale fra i due segnali:
% % 
% % Nella parte iniziale (in basso a sinistra della matrice):
% % active → anticipa → obs (correlazioni positive alte).
% % 
% % Nella parte finale (in alto a destra):
% % active → termina prima di → obs, che rimane alto più a lungo (correlazioni negative alte).
% 
% 
% 
% %%%%%%%%%%%%%%%%% punti 3-4-5 nel paragrafo Model-based brain-to-brain coupling
% % MEtodi (pag e4) we develop a framework for assessing five types of generalization simultaneously:
% % testing encoding model generalization
% % 3)across different brain regions (e.g., from SM to STG electrodes), 
% % 4) across tasks/processes (speaking/production and listening/comprehension)
% % 5) across lags (e.g., speaker pre-word onset to listener post-word onset).
% 
% %% sostanzialmente fa la stessa cosa (calcolo correlazioni, media delle stesse etc)
% % cambiando però i riferiemnti; nel punto 3 fa la stessa cosa ma usando
% % aree diverse del cervello (stima quindi il modello considerando alcuni
% % elettrodi e vede la capacità predittiva su altri elettrodi
% 
% 
