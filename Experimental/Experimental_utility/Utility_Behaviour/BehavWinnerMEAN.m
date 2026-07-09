function [Data_result,TableSolo,TableJoint] = BehavWinnerMEAN(dataSession, par)

mode = par.mode;
num_dir = par.num_dir;

Data_result = struct();
VarNames = {'Direction','Winner','NormalityTest','Hnorm','testType','pvalue','Significant','S_mean','K_mean'};
VarTypes = {'double','string','string','double','string','double','string','double','double'};
TableSolo = table('Size',[num_dir,length(VarNames)],...
    'VariableNames', VarNames,'VariableTypes',VarTypes);
TableJoint = table('Size',[num_dir,length(VarNames)],...
    'VariableNames', VarNames,'VariableTypes',VarTypes);

for dirIdx = 1:num_dir
    dirName = sprintf('dir%d', dirIdx);
    S1_mean = NaN(length(dataSession),1);
    K2_mean = NaN(length(dataSession),1);
    JS3_mean = NaN(length(dataSession),1);
    JK3_mean = NaN(length(dataSession),1);

    for nsess = 1:length(dataSession)
        data_sess = dataSession(nsess).cond;

        % Extraction
        S1_mean(nsess,1) = data_sess(1).(dirName); % Solo S
        K2_mean(nsess,1) = data_sess(2).(dirName); % Solo K
        JS3_mean(nsess,1) = data_sess(3).(dirName); % Joint S
        JK3_mean(nsess,1) = data_sess(4).(dirName); % Joint K
    end

    if strcmp(mode,'time')
        if mean(S1_mean) < mean(K2_mean)
            winnerSolo = "S";
        elseif mean(S1_mean) > mean(K2_mean)
            winnerSolo = "K";
        else
            winnerSolo = "tie";
        end

        if mean(JS3_mean) < mean(JK3_mean)
            winnerJoint = "S";
        elseif mean(JS3_mean) > mean(JK3_mean)
            winnerJoint = "K";
        else
            winnerJoint = "tie";
        end
    else
        if mean(S1_mean) > mean(K2_mean)
            winnerSolo = "S";
        elseif mean(S1_mean) < mean(K2_mean)
            winnerSolo = "K";
        else
            winnerSolo = "tie";
        end

        if mean(JS3_mean) > mean(JK3_mean)
            winnerJoint = "S";
        elseif mean(JS3_mean) < mean(JK3_mean)
            winnerJoint = "K";
        else
            winnerJoint = "tie";
        end
    end
    meanS1 = mean(S1_mean);
    meanK2 = mean(K2_mean);
    meanJS3 = mean(JS3_mean);
    meanJK3 = mean(JK3_mean);

    Data_result.Solo.(dirName)      = winnerSolo;
    Data_result.Joint.(dirName)     = winnerJoint;

    %% Normality and statistic
    dSolo = S1_mean - K2_mean;
    n = length(dSolo);
    if n < 50
        % Shapiro-Wilk
        [HnormSolo, pNormSolo] = swtest(dSolo);
        NormalityTestSolo = 'Shapiro–Wilk';
    else
        % Kolmogorov-Smirnov
        z = (dSolo - mean(dSolo)) / std(dSolo);
        [HnormSolo, pNormSolo] = kstest(z);
        NormalityTestSolo = 'Kolmogorov–Smirnov';
    end

    if HnormSolo == 0
        [hSolo, pValSolo] = ttest(S1_mean, K2_mean);
        testTypeSolo = 'paired t-test';
    else
        % Non normale → Wilcoxon signed-rank
        pValSolo = signrank(S1_mean, K2_mean);
        testTypeSolo = 'Wilcoxon signed-rank';
    end

    if pValSolo < 0.001
        signifSolo = '***';
    elseif pValSolo < 0.01
        signifSolo = '**';
    elseif pValSolo < 0.05
        signifSolo = '*';
    else
        signifSolo = 'None';
    end

    TableSolo.Direction(dirIdx)     = dirIdx;
    TableSolo.Winner(dirIdx)        = winnerSolo;
    TableSolo.NormalityTest(dirIdx) = NormalityTestSolo;
    TableSolo.Hnorm(dirIdx)         = HnormSolo;
    TableSolo.testType(dirIdx)      = testTypeSolo;
    TableSolo.pvalue(dirIdx)        = pValSolo;
    TableSolo.Significant(dirIdx)   = signifSolo;
    TableSolo.S_mean(dirIdx)        = meanS1;
    TableSolo.K_mean(dirIdx)        = meanK2;

    dJoint = JS3_mean - JK3_mean;
    n = length(dJoint);
    if n < 50
        % Shapiro-Wilk
        [HnormJoint, pNormJoint] = swtest(dJoint);
        NormalityTestJoint = 'Shapiro–Wilk';
    else
        % Kolmogorov-Smirnov
        z = (dJoint - mean(dJoint)) / std(dJoint);
        [HnormJoint, pNormJoint] = kstest(z);
        NormalityTestJoint = 'Kolmogorov–Smirnov';
    end

    if HnormJoint == 0
        [hJoint, pValJoint] = ttest(JS3_mean, JK3_mean);
        testTypeJoint = 'paired t-test';
    else
        % Non normale → Wilcoxon signed-rank
        pValJoint = signrank(JS3_mean, JK3_mean);
        testTypeJoint = 'Wilcoxon signed-rank';
    end

    if pValJoint < 0.001
        signifJoint = '***';
    elseif pValJoint < 0.01
        signifJoint = '**';
    elseif pValJoint < 0.05
        signifJoint = '*';
    else
        signifJoint = 'None';
    end

    TableJoint.Direction(dirIdx)     = dirIdx;
    TableJoint.Winner(dirIdx)        = winnerJoint;
    TableJoint.NormalityTest(dirIdx) = NormalityTestJoint;
    TableJoint.Hnorm(dirIdx)         = HnormJoint;
    TableJoint.testType(dirIdx)      = testTypeJoint;
    TableJoint.pvalue(dirIdx)        = pValJoint;
    TableJoint.Significant(dirIdx)   = signifJoint;
    TableJoint.S_mean(dirIdx)        = meanJS3;
    TableJoint.K_mean(dirIdx)        = meanJK3;

end