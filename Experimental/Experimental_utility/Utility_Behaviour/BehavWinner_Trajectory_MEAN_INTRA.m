function [Data_result, TableSoloSvsJointS, TableSoloKvsJointK] = BehavWinner_Trajectory_MEAN_INTRA(data_S, data_K, par)
% Funzione DUALE: Statistica SoloS vs JointS e SoloK vs JointK
% Confronti INTRA-monkey (S con S, K con K)
num_dir = par.num_dir;
inField = par.InField; % 'angular_error', 'exitTime', 'ExitPerf'
num_sessions = length(data_S);
Data_result = struct();

% --- SETUP TABELLE ---
VarNames = {'Metric', 'Direction','Winner','NormalityTest','Hnorm','testType','pvalue','Significant','Solo_mean','Joint_mean'};
VarTypes = {'string', 'double','string','string','double','string','double','string','double','double'};

TableSoloSvsJointS = table('Size',[num_dir,length(VarNames)], 'VariableNames', VarNames, 'VariableTypes', VarTypes);
TableSoloKvsJointK = table('Size',[num_dir,length(VarNames)], 'VariableNames', VarNames, 'VariableTypes', VarTypes);

for dirIdx = 1:num_dir
    dirName = sprintf('dir%d', dirIdx);
    
    % Preallocazione vettori medie per sessione
    S1vsJS3_sess_means = NaN(num_sessions, 2);  % col1: SoloS, col2: JointS
    K2vsJK3_sess_means = NaN(num_sessions, 2);  % col1: SoloK, col2: JointK
    
    for nsess = 1:num_sessions
        % SESSIONE S: SoloS(1) vs JointS(3)
        vS1_raw = data_S(nsess).(inField)(1).(dirName);
        vJS3_raw = data_S(nsess).(inField)(3).(dirName);
        
        % SESSIONE K: SoloK(2) vs JointK(3)
        vK2_raw = data_K(nsess).(inField)(2).(dirName);
        vJK3_raw = data_K(nsess).(inField)(3).(dirName);

        S1vsJS3_sess_means(nsess,1) = mean(vS1_raw);
        S1vsJS3_sess_means(nsess,2) = mean(vJS3_raw);
        K2vsJK3_sess_means(nsess,1) = mean(vK2_raw);
        K2vsJK3_sess_means(nsess,2) = mean(vJK3_raw);

    end

    % --- S: SoloS vs JointS ---
    mS1  = mean(S1vsJS3_sess_means(:,1));
    mJS3 = mean(S1vsJS3_sess_means(:,2));

    if strcmp(inField, 'ExitPerf') % Più alto meglio
        if mS1 > mJS3
            winnerS = "SoloS";
        else
            winnerS = "JointS";
        end
    else % Più basso meglio
        if mS1 < mJS3
            winnerS = "SoloS";
        else
            winnerS = "JointS";
        end
    end
    
    [~, pS, testS, normS, hNormS] = performComparison(...
        S1vsJS3_sess_means(:,1), S1vsJS3_sess_means(:,2));
    sigS = getSigStars(pS);
    
    TableSoloSvsJointS(dirIdx,:) = {inField, dirIdx, winnerS, normS, hNormS, testS, pS, sigS, mS1, mJS3};
    
    % --- K: SoloK vs JointK ---
    mK2  = mean(K2vsJK3_sess_means(:,1));
    mJK3 = mean(K2vsJK3_sess_means(:,2));
    
    if strcmp(inField, 'ExitPerf')
        if mK2 > mJK3
            winnerK = "SoloK";
        else
            winnerK = "JointK";
        end
    else
        if mK2 < mJK3
            winnerK = "SoloK";
        else
            winnerK = "JointK";
        end
    end
    
    [~, pK, testK, normK, hNormK] = performComparison(...
        K2vsJK3_sess_means(:,1), K2vsJK3_sess_means(:,2));
    sigK = getSigStars(pK);
    
    TableSoloKvsJointK(dirIdx,:) = {inField, dirIdx, winnerK, normK, hNormK, testK, pK, sigK, mK2, mJK3};
    
    % Output Struct
    Data_result.S.(dirName).Winner  = winnerS;
    Data_result.S.(dirName).pValue  = pS;
    Data_result.K.(dirName).Winner  = winnerK;
    Data_result.K.(dirName).pValue  = pK;
end
end

%% FUNZIONI DI SUPPORTO (invariate)
function [h, p, testType, normType, hNorm] = performComparison(vec1, vec2)
    diffVal = vec1 - vec2;
    n = length(diffVal);
    
    if n < 50
        [hNorm, ~] = swtest(diffVal);
        normType = 'Shapiro–Wilk';
    else
        z = (diffVal - mean(diffVal)) / std(diffVal);
        [hNorm, ~] = kstest(z);
        normType = 'Kolmogorov–Smirnov';
    end
    
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
