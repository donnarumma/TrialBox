function [Data_result, TableSoloSvsJointS, TableSoloKvsJointK] = BehavWinnerMEAN_INTRA(dataSession, par)
% DUALE: SoloS vs JointS e SoloK vs JointK (INTRA-monkey)
mode = par.mode;
num_dir = par.num_dir;
Data_result = struct();

VarNames = {'Direction','Winner','NormalityTest','Hnorm','testType','pvalue','Significant','Solo_mean','Joint_mean'};
VarTypes = {'double','string','string','double','string','double','string','double','double'};

TableSoloSvsJointS = table('Size',[num_dir,length(VarNames)],...
    'VariableNames', VarNames,'VariableTypes',VarTypes);
TableSoloKvsJointK = table('Size',[num_dir,length(VarNames)],...
    'VariableNames', VarNames,'VariableTypes',VarTypes);

for dirIdx = 1:num_dir
    dirName = sprintf('dir%d', dirIdx);
    
    % S: SoloS(1) vs JointS(3)
    S1vsJS3_mean = NaN(length(dataSession),2);
    % K: SoloK(2) vs JointK(4)
    K2vsJK4_mean = NaN(length(dataSession),2);
    
    for nsess = 1:length(dataSession)
        data_sess = dataSession(nsess).cond;
        
        % SoloS vs JointS
        S1vsJS3_mean(nsess,1) = data_sess(1).(dirName); % SoloS
        S1vsJS3_mean(nsess,2) = data_sess(3).(dirName); % JointS
        
        % SoloK vs JointK
        K2vsJK4_mean(nsess,1) = data_sess(2).(dirName); % SoloK
        K2vsJK4_mean(nsess,2) = data_sess(4).(dirName); % JointK
    end
    
    %% S: SoloS vs JointS
    if strcmp(mode,'time') % Più basso meglio
        if mean(S1vsJS3_mean(:,1)) < mean(S1vsJS3_mean(:,2))
            winnerS = "SoloS";
        elseif mean(S1vsJS3_mean(:,1)) > mean(S1vsJS3_mean(:,2))
            winnerS = "JointS";
        else
            winnerS = "tie";
        end
    else % Più alto meglio
        if mean(S1vsJS3_mean(:,1)) > mean(S1vsJS3_mean(:,2))
            winnerS = "SoloS";
        elseif mean(S1vsJS3_mean(:,1)) < mean(S1vsJS3_mean(:,2))
            winnerS = "JointS";
        else
            winnerS = "tie";
        end
    end
    
    %% K: SoloK vs JointK
    if strcmp(mode,'time')
        if mean(K2vsJK4_mean(:,1)) < mean(K2vsJK4_mean(:,2))
            winnerK = "SoloK";
        elseif mean(K2vsJK4_mean(:,1)) > mean(K2vsJK4_mean(:,2))
            winnerK = "JointK";
        else
            winnerK = "tie";
        end
    else
        if mean(K2vsJK4_mean(:,1)) > mean(K2vsJK4_mean(:,2))
            winnerK = "SoloK";
        elseif mean(K2vsJK4_mean(:,1)) < mean(K2vsJK4_mean(:,2))
            winnerK = "JointK";
        else
            winnerK = "tie";
        end
    end
    
    % Medie finali
    meanS1_SoloS  = mean(S1vsJS3_mean(:,1));
    meanJS3_JointS= mean(S1vsJS3_mean(:,2));
    meanK2_SoloK  = mean(K2vsJK4_mean(:,1));
    meanJK4_JointK= mean(K2vsJK4_mean(:,2));
    
    Data_result.S.(dirName)   = winnerS;
    Data_result.K.(dirName)   = winnerK;
    
    %% STATISTICA SoloS vs JointS
    dS = S1vsJS3_mean(:,1) - S1vsJS3_mean(:,2);
    [HnormS, ~] = normalityTest(dS);
    [pValS, testTypeS, signifS] = statTest(S1vsJS3_mean(:,1), S1vsJS3_mean(:,2), HnormS);
    
    TableSoloSvsJointS(dirIdx,:) = {dirIdx, winnerS, HnormS.testType, HnormS.h, testTypeS, pValS, signifS, meanS1_SoloS, meanJS3_JointS};
    
    %% STATISTICA SoloK vs JointK
    dK = K2vsJK4_mean(:,1) - K2vsJK4_mean(:,2);
    [HnormK, ~] = normalityTest(dK);
    [pValK, testTypeK, signifK] = statTest(K2vsJK4_mean(:,1), K2vsJK4_mean(:,2), HnormK);
    
    TableSoloKvsJointK(dirIdx,:) = {dirIdx, winnerK, HnormK.testType, HnormK.h, testTypeK, pValK, signifK, meanK2_SoloK, meanJK4_JointK};
end
end

%% HELPER FUNCTIONS (riutilizzabili)
function [Hnorm, pNorm] = normalityTest(data)
    n = length(data);
    if n < 50
        [Hnorm.h, pNorm] = swtest(data);
        Hnorm.testType = 'Shapiro–Wilk';
    else
        z = (data - mean(data)) / std(data);
        [Hnorm.h, pNorm] = kstest(z);
        Hnorm.testType = 'Kolmogorov–Smirnov';
    end
end

function [pVal, testType, signif] = statTest(vec1, vec2, Hnorm)
    if Hnorm.h == 0
        [~, pVal] = ttest(vec1, vec2);
        testType = 'paired t-test';
    else
        pVal = signrank(vec1, vec2);
        testType = 'Wilcoxon signed-rank';
    end
    signif = getSigStars(pVal);
end

function stars = getSigStars(p)
    if p < 0.001, stars = '***';
    elseif p < 0.01, stars = '**';
    elseif p < 0.05, stars = '*';
    else, stars = 'None';
    end
end
