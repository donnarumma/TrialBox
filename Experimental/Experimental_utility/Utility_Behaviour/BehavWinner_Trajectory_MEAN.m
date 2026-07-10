function [Data_result, TableSolo, TableJoint] = BehavWinner_Trajectory_MEAN(data_S, data_K, par)
% Funzione 4: Statistica sulle medie per dati Jxy (Traiettorie)
% Duale alla Funzione 3, ma strutturata per l'input Jxy della Funzione 2.

num_dir = par.num_dir;
inField = par.InField; % es. 'angular_error', 'exitTime', 'ExitPerf'
num_sessions = length(data_S);

Data_result = struct();

% --- SETUP TABELLE (Stile Funzione 3) ---
VarNames = {'Metric', 'Direction','Winner','NormalityTest','Hnorm','testType','pvalue','Significant','S_mean','K_mean'};
VarTypes = {'string', 'double','string','string','double','string','double','string','double','double'};

TableSolo = table('Size',[num_dir,length(VarNames)], 'VariableNames', VarNames, 'VariableTypes', VarTypes);
TableJoint = table('Size',[num_dir,length(VarNames)], 'VariableNames', VarNames, 'VariableTypes', VarTypes);

for dirIdx = 1:num_dir
    dirName = sprintf('dir%d', dirIdx);
    
    % Preallocazione vettori per le medie di ogni sessione
    S1_sess_means = NaN(num_sessions, 1);
    K2_sess_means = NaN(num_sessions, 1);
    JS3_sess_means = NaN(num_sessions, 1);
    JK3_sess_means = NaN(num_sessions, 1);
    
    for nsess = 1:num_sessions
        % Estrazione dati raw per la sessione corrente
        vS1_raw = data_S(nsess).(inField)(1).(dirName);
        vK2_raw = data_K(nsess).(inField)(2).(dirName);
        vJS3_raw = data_S(nsess).(inField)(3).(dirName);
        vJK3_raw = data_K(nsess).(inField)(3).(dirName);
        
        % Calcolo media sessione
        S1_sess_means(nsess)  = mean(vS1_raw);
        K2_sess_means(nsess)  = mean(vK2_raw);
        JS3_sess_means(nsess) = mean(vJS3_raw);
        JK3_sess_means(nsess) = mean(vJK3_raw);

    end
    
    % --- DETERMINAZIONE VINCITORE (Logica Funzione 2) ---
    mS1 = mean(S1_sess_means);  mK2 = mean(K2_sess_means);
    mJS3 = mean(JS3_sess_means); mJK3 = mean(JK3_sess_means);
    
    if strcmp(inField, 'ExitPerf') % Più alto è meglio
        if mS1 > mK2, winnerSolo = "S"; else, winnerSolo = "K"; end
        if mJS3 > mJK3, winnerJoint = "S"; else, winnerJoint = "K"; end
    else % Più basso è meglio (Errori, Tempi)
        if mS1 < mK2, winnerSolo = "S"; else, winnerSolo = "K"; end
        if mJS3 < mJK3, winnerJoint = "S"; else, winnerJoint = "K"; end
    end
    
    % --- ANALISI STATISTICA (Logica Funzione 3) ---
    
    % 1. Analisi SOLO
    [~, pS, testS, normS, hNormS] = performComparison(S1_sess_means, K2_sess_means);
    sigS = getSigStars(pS);
    
    TableSolo(dirIdx,:) = {inField, dirIdx, winnerSolo, normS, hNormS, testS, pS, sigS, mS1, mK2};
    
    % 2. Analisi JOINT
    [~, pJ, testJ, normJ, hNormJ] = performComparison(JS3_sess_means, JK3_sess_means);
    sigJ = getSigStars(pJ);
    
    TableJoint(dirIdx,:) = {inField, dirIdx, winnerJoint, normJ, hNormJ, testJ, pJ, sigJ, mJS3, mJK3};
    
    % Output Struct
    Data_result.Solo.(dirName).Winner = winnerSolo;
    Data_result.Solo.(dirName).pValue = pS;
    Data_result.Joint.(dirName).Winner = winnerJoint;
    Data_result.Joint.(dirName).pValue = pJ;
end

end

%% --- FUNZIONI DI SUPPORTO (Helpers) ---

function [h, p, testType, normType, hNorm] = performComparison(vec1, vec2)
    diffVal = vec1 - vec2;
    n = length(diffVal);
    
    % Test Normalità
    if n < 50
        [hNorm, ~] = swtest(diffVal);
        normType = 'Shapiro–Wilk';
    else
        z = (diffVal - mean(diffVal)) / std(diffVal);
        [hNorm, ~] = kstest(z);
        normType = 'Kolmogorov–Smirnov';
    end
    
    % Scelta Test Statistico
    if hNorm == 0
        [h, p] = ttest(vec1, vec2);
        testType = 'paired t-test';
    else
        p = signrank(vec1, vec2);
        h = p <= 0.05;
        testType = 'Wilcoxon signed-rank';
    end
end

function stars = getSigStars(p)
    if p < 0.001, stars = '***';
    elseif p < 0.01, stars = '**';
    elseif p < 0.05, stars = '*';
    else, stars = 'None';
    end
end