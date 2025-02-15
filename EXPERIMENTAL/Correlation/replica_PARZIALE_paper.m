%%% Ref libreria python https://zenodo.org/records/12359299

% Replica moooolto semplificata del paper
% mi sono focalizzato sui punti chiave; 
% 1) ho simulato i dati (quindi saltato il signal processing e
% l'estrazione di embedding con gpt2);
% 2) ho allineato i dati campionando come lo fanno loro
% 3) regressione ridge per lag e per elettrodo con cross validation sui dati 
%    e calcolo correlazioni tra valori predetti e osservati - (correlazione
%    tra segmenti dello stimolo)
% 4) (Salto la parte di selezione degli elettrodi) Brain to brain coupling 
%  praticamente testano il modello allenato su un soggetto, sull'altro
%  soggetto...

%%% Le parti di riferimento sul paper sono nei metodi del paper:
% "A shared model-based linguistic space for transmitting our thoughts from brain to brain in natural conversations"
% Encoding analysis
% Model-based brain-to-brain coupling


%%%%%%%%%%%%%%%%%%%%%%%%%% GENERAZIONE DATI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generiamo dei simil-dati dialogo per due soggetti (durata un minuto)
%% neurali 
% frequenza campionamento in hertz
fs=512;
% durata dialogo 60 secondi
dd=60;
% tempo 
tt = linspace(0,dd, dd*fs);
% elettrodi
n_channels=64;
% parole al minuto
wpm=150;

%% genero dati EcoG per soggetto (elettrodi per tempo con componente oscillatoria
%% a 40 hert e rumore gaussiano)
% speaker
speaker_neur= 0.5 * sin(2 * pi * 40 * tt)' * ones(1, n_channels) + randn(length(tt), n_channels);
%listener
listener_neur=sin(2 * pi * 30 * tt)' * ones(1, n_channels) + sin(2 * pi * 60 * tt)' * ones(1, n_channels) + ...
       randn(length(tt), n_channels) * 0.2;

% Onset casuali delle parole distribuiti nel minuto di dialogo (in
% frequenza...evito i primi due e gli ultimi due secondi per sampling successivo)
% ampiezza intorno degli onset delle parole
intorno_sec=4;
word_onsets = sort(randi([intorno_sec/2*fs length(tt)-intorno_sec/2*fs], 1, wpm));

% embeddings - casuali - tipo GPT2 (che per ogni parola tira fuori un embeddding di
% 1600 componenti)...qui metto meno componenti
embedding_size = 500; 
embeddings = randn(wpm, embedding_size); 


%% %%%%%%&&%%%%%%%% ALLINEAMENTO DATI NEURALI - DATI PAROLE %%%%%%%%%%%%%%%%
% di fatto si analizza l'attività neurale in intervalli di tempo specifici
% intorno agli onsets delle parole (loro lo fanno 4 secondi prima e 4 dopo
% ogni parola, qui lo faccio di due)

% ogni intorno di 4 secondi viene binnato in finestre di dimensione 
% arbitraria (250ms' con overlap 
% di 62.5ms ...quindi se una parola inizia a 5.5 sec. il suo intorno sarà tra 
% 3.5 e 7.5 sec)
% Definizione dei lags espressi in numero di campioni
ww=250; %  finestra in ms
% jump per ogni blocco
overlap=62.5  %ms
% % Conversione della finestra e dell'overlap da millisecondi a numero di campioni
n_samples = (fs/ 1000*ww);
overlap_samples = (fs/ 1000*overlap);
% Numero totale di blocchi nel periodo di interesse
num_blocchi = (intorno_sec * fs) / (ww);
 
for i = 1:wpm
    idx_min = word_onsets(i) - (intorno_sec/2 * fs);
    idx_max = word_onsets(i) + (intorno_sec/2 * fs);

    start_idx = idx_min;
    blocchi_idx = []; 
    
    while start_idx + n_samples - 1 <= idx_max
        end_idx = start_idx + n_samples - 1;
        blocchi_idx = [blocchi_idx; start_idx, end_idx];
        
        % Calcola il prossimo indice di partenza con l'overlap
        start_idx = start_idx + overlap_samples;
    end

    %  alla fine mi ritrovo con due strutture (per speaker e listener)
    % che rappresentano l'attività neurale dei soggetti realtiva alle
    % parole
    % per cui avremo un tot (61) di  blocchi per ogni parola
    % di dimensione 128x64 (obs in hz * elettrodi)
    dati_speaker{i} = cell(size(blocchi_idx, 1), 1);
    dati_listener{i} = cell(size(blocchi_idx, 1), 1);
    for j = 1:size(blocchi_idx, 1)
        start_blocco = blocchi_idx(j, 1);
        end_blocco = blocchi_idx(j, 2);
        dati_speaker{i}{j} = speaker_neur(start_blocco:end_blocco, :);
        dati_listener{i}{j} =listener_neur(start_blocco:end_blocco, :);

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% RIDGE REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% la ridge regression viene fatta per ogni blocco e per ogni elettrodo
% devo quindi estrarre dalla struttura i vettori (128*1) per parola per
% elettrodo (nota che dal punto precdente ricaviamo che ogni 

% Numero di blocchi per parola
n_blocchi = size(dati_speaker{1}, 1); 
n_channels % Numero di elettrodi per blocco
n_samples % lunghezza blocco

% struttura Y (speaker e listener)
% (saranno tanti vettori di dimensione (numero parole *1) 
% un vettore per ogni lag per ogni elettrodo...siccome ogni lag (per
% ogni elettrodo) include un blocco di osservazioni per parola, si fa una
% media sul blocco...quidni un valore per ogni parola per ogni lag per ogni elettrodo
% Blocchi di variabili di rispsota per soggetto

Y_s= cell(n_blocchi, n_channels);
Y_l= cell(n_blocchi, n_channels);
% Riempimento di Y con i dati
for b = 1:n_blocchi
    for e = 1:n_channels
        Y_s{b,e} = [];
        Y_l{b,e} = [];
        for p = 1:wpm
            
            Y_s{b,e} = [Y_s{b,e}; mean(dati_speaker{p}{b}(:, e))];
            Y_l{b,e} = [Y_l{b,e}; mean(dati_listener{p}{b}(:, e))];
        end
    end
end

% quindi ora abbiamo due strutture Y che contengono  61*64 elementi
% (blocchi (o lags)*elettrodi) da 150 parole
% osservazioni (medie per lag per parola). 
% recuperiamo la matrce X (embeddig per parola)

X=embeddings
 
%% Regressioni ridge 
lambda = 1000; % Supponiamo un lambda già ottimale

% Inizializzazione delle strutture per i coefficienti e predizioni
% speaker
coefs_s = cell(n_blocchi, n_channels);
y_hat_s = cell(n_blocchi, n_channels);
% listener
coefs_l = cell(n_blocchi, n_channels);
y_hat_l = cell(n_blocchi, n_channels);

%%% ottengo quindi delle y stimate e dei coefficienti (uno per ogni aspetto
%%% dell'embedding delle parole...quindi avendo generato  embedding di 500
%%% valori avrò 500 coefficienti per regressione)
for b = 1:n_blocchi
    for e = 1:n_channels
        % target (segnale neurale per quel blocco ed elettrodo)
        y_s = Y_s{b, e}; 
        y_l = Y_l{b, e};
        %%% check su dimensioni
        %if size(X,1) ~= size(Y_target,1)
         %   error('Dimensioni non corrispondenti tra X e y target in blocco %d, elettrodo %d', b, e);
        %end

        % Stima dei coefficienti della Ridge Regression: (X'X + lambda*I) \
        % X'y
        % Matrice identità per la regolarizzazione
        I = eye(size(X,2)); 
        coefs_s{b, e} = (X' * X + lambda * I) \ (X' * y_s);
        coefs_l{b, e} = (X' * X + lambda * I) \ (X' * y_l);
        % Predizione della risposta neurale
        y_hat_s{b, e} = X * coefs_s{b, e};
        y_hat_l{b, e} = X * coefs_l{b, e};

    end
end

%% %%%%% Correlation
% Questo di fatto è il punto 1 nel paragrafo Model-based brain-to-brain
% coupling  Metodi (pag e4)
% "we develop a framework for assessing five types of generalization simultaneously:
% testing encoding model generalization
% (1) across segments of the stimulus (using 10-fold cross-validation)"
%
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


%% faccio 3 cv
% Inizializzazione delle strutture per i coefficienti, predizioni e
% correlazioni  
K           = 3; % Numero di fold per la cross-validation 
coefs_s_cv  = cell(n_blocchi, n_channels);
y_hat_s_cv  = cell(n_blocchi, n_channels);
corrs_s_cv  = zeros(n_blocchi, n_channels, K);
%listener
coefs_l_cv  = cell(n_blocchi, n_channels);
y_hat_l_cv  = cell(n_blocchi, n_channels);
corrs_l_cv  = zeros(n_blocchi, n_channels, K);
foldIndices = ceil(linspace(1, wpm+1, K+1)); 


for b = 1:n_blocchi
    for e = 1:n_channels
        % y osservate speaker
        % Segnale neurale target per il blocco b, elettrodo e
        y_s = Y_s{b, e};
        % y osservate listener
        y_l = Y_l{b, e};
        for k = 1:K
            % Definisco train/test
            testIdx = foldIndices(k):foldIndices(k+1)-1;
            trainIdx = setdiff(1:wpm, testIdx);
            
            X_train = X(trainIdx, :);
            y_train_s = y_s(trainIdx, :);
            y_train_l = y_l(trainIdx, :);

            X_test = X(testIdx, :);
            y_test_s = y_s(testIdx, :);
            y_test_l = y_l(testIdx, :);

           
            % Stima dei coefficienti della Ridge Regression: (X'X + lambda*I) \
            % X'y
            % Matrice identità per la regolarizzazione
            I = eye(size(X,2));
            % speaker
            coefs_s_cv{b, e, k} = (X_train' * X_train + lambda * I) \ (X_train' * y_train_s);
            %listener
            coefs_l_cv{b, e, k} = (X_train' * X_train + lambda * I) \ (X_train' * y_train_l);

            % Predicted y
            y_hat_test_s_cv = X_test * coefs_s_cv{b, e, k};
            y_hat_test_l_cv = X_test * coefs_l_cv{b, e, k};

            
            % correlazione predicted observed (speaker)
            corrs_s_cv(b, e, k) = corr(y_test_s, y_hat_test_s_cv);
            corrs_l_cv(b, e, k) = corr(y_test_l, y_hat_test_s_cv);
        end
    end
end

% Media della correlazione speaker
% alla fine mi ritrovo con una matrice di correlazione per blocco per 
% per elettrodo
corrs_mean_s_cv = mean(corrs_s_cv, 3);
% listener
corrs_mean_l_cv = mean(corrs_l_cv, 3);


% Display risultati
disp("Correlazione media per ogni canale e blocco - within subject - speaker :");
disp(corrs_mean_s_cv);
disp("Correlazione media per ogni canale e blocco - within subject - listener :");
disp(corrs_mean_l_cv);
%%%%%%%%%%%%%%%%%%%%%%% 4) Brain  to Brain coupling %%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%% punto 2 nel paragrafo Model-based brain-to-brain coupling MEtodi (pag e4)
%%%%%%%%%%%%%%%%%%%%  GENERALIZATION ACROSS SUBJECTS ######################
% "we develop a framework for assessing five types of generalization simultaneously:
% testing encoding model generalization
% (2) across subjects (within speaker–listener dyads)"
% considero i coefficienti stimati  dello speaker per capire se l'attività
% neurale dello speaker è predittiva per  il listener e viceversa e se la
% relazione tra le due attività cerebrali segue schemi coerenti 

% di fatto le y previste non cambiano da quelle stimate nel modello within-
% subject...si usano le stesse X (embedding parole) e stessi coefficienti

% 3 fold cv
%% 

K=3
%Inizializzazione della struttura per la correlazione inter-soggetto
% speaker-listener
corrs_inter_s_l = zeros(n_blocchi, n_channels, K);
% listener-speaker
corrs_inter_l_s = zeros(n_blocchi, n_channels, K);
for b = 1:n_blocchi
    for e = 1:n_channels
        y_listener = Y_l{b, e}; % Dati neurali del listener
        y_speaker = Y_s{b, e}; % Dati neurali del listener
        
        for k = 1:K
            % Indici del test (come sopra)
            testIdx = foldIndices(k):foldIndices(k+1)-1;
            
            % Dati test del listener
            X_test = X(testIdx, :);
            y_test_listener = y_listener(testIdx, :);
            y_test_speaker  = y_speaker(testIdx, :);

            % Uso i coefficienti dello speaker sul listener e viceversa
            % senza riadattare!
            y_hat_listener  = X_test*coefs_s_cv{b, e, k};
            y_hat_speaker   = X_test*coefs_l_cv{b, e, k};

            % Correlazione tra osservato (listener) e predetto (da speaker)
            corrs_inter_s_l(b, e, k) = corr(y_test_listener, y_hat_listener);
            corrs_inter_l_s(b, e, k) = corr(y_test_speaker, y_hat_speaker);

        end
    end
end

% di fatto mi trovo con delle  matrici 3d con le correlazioni per elettrodi 
% per lags per folds [n_blocchi, n_channels, n_folds]

%% %%%% COSTRUZIONE HEATMAPS lag lag fig 3a %%%%%%%%%%
n_lags=n_blocchi; n_folds=K; 

% Calcolo la media across electrodes per ogni fold (come fanno in Python)
%  ottengo matrici (lags, folds)
corrs_s_l_avg = squeeze(mean(corrs_inter_s_l, 2)); 
corrs_l_s_avg = squeeze(mean(corrs_inter_l_s, 2));

% Offset per evitare problemi di numeri molto piccoli
epsilon = 1e-4; % Valore piccolo per evitare NaN nelle correlazioni

% Aggiungiamo epsilon ai valori troppo vicini a zero
corrs_s_l_avg(abs(corrs_s_l_avg) < epsilon) = sign(corrs_s_l_avg(abs(corrs_s_l_avg) < epsilon)) * epsilon;
corrs_l_s_avg(abs(corrs_l_s_avg) < epsilon) = sign(corrs_l_s_avg(abs(corrs_l_s_avg) < epsilon)) * epsilon;

% Inizializzo la matrice lag-lag
corr_matrix = zeros(n_lags, n_lags);

% Calcolo  l correlazione per ogni combinazione di lag speaker-listener
for i = 1:n_lags
    for j = 1:n_lags
        % Vettore di 3 valori per il lag i (Speaker → Listener)
        speaker_vector = corrs_s_l_avg(i, :); % (1 × 3)
        
        % Vettore di 3 valori per il lag j (Listener → Speaker)
        listener_vector = corrs_l_s_avg(j, :); % (1 × 3)
        
        % Calcolo la correlazione tra i due vettori
        corr_matrix(i, j) = corr(speaker_vector', listener_vector');
    end
end


% ora medio la matrice di correlazione across folds
% corr_matrix = mean(corr_matrix_folds, 3);

%% per il plot
% time range
time_range = linspace(-2, 2, n_lags); % Array dei tempi associati ai lag
tick_indices = round(linspace(1, n_lags, 5)); % Seleziona 5 tick equidistanti
tick_labels = arrayfun(@(x) sprintf('%.1fs', time_range(x)), tick_indices, 'UniformOutput', false);

% Visualizzazione della heatmap
figure;
imagesc(corr_matrix);
colormap('hot');
colorbar;
axis square;
xlabel('Lag Listener (s)');
ylabel('Lag Speaker (s)');
title('Heatmap della correlazione lag-lag tra speaker e listener');
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_labels, ...
         'YTick', tick_indices, 'YTickLabel', tick_labels);

%%%%%%%%%%%%%%%%% punti 3-4-5 nel paragrafo Model-based brain-to-brain coupling
% MEtodi (pag e4) we develop a framework for assessing five types of generalization simultaneously:
% testing encoding model generalization
% 3)across different brain regions (e.g., from SM to STG electrodes), 
% 4) across tasks/processes (speaking/production and listening/comprehension)
% 5) across lags (e.g., speaker pre-word onset to listener post-word onset).

%% sostanzialmente fa la stessa cosa (calcolo correlazioni, media delle stesse etc)
% cambiando però i riferiemnti; nel punto 3 fa la stessa cosa ma usando
% aree diverse del cervello (stima quindi il modello considerando alcuni
% elettrodi e vede la capacità predittiva su altri elettrodi


