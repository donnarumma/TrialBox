function [Data_result,BinomialTableSolo,BinomialTableJoint] = BehavWinner(data, par)

mode = par.mode;
num_dir = par.num_dir;

Data_result = struct();
VarNames = {'Direction','Wins_S','Wins_K','Winner','p_binom','Significant','S_Percent','K_Percent'};
VarTypes = {'double','double','double','string','double','double','double','double'};
BinomialTableSolo = table('Size',[num_dir,length(VarNames)],...
    'VariableNames', VarNames,'VariableTypes',VarTypes);
BinomialTableJoint = table('Size',[num_dir,length(VarNames)],...
    'VariableNames', VarNames,'VariableTypes',VarTypes);
for dirIdx = 1:num_dir
    dirName = sprintf('dir%d', dirIdx);
    wins_S1 = 0;
    wins_K2 = 0;

    wins_JS3 = 0;
    wins_JK3 = 0;
    for nsess = 1:length(data)
        data_sess = data(nsess).cond;

        % Extraction
        S1 = data_sess(1).(dirName); % Solo S
        K2 = data_sess(2).(dirName); % Solo K
        JS3 = data_sess(3).(dirName); % Joint S
        JK3 = data_sess(4).(dirName); % Joint K

        % mean
        S1_mean     = mean(S1);
        K2_mean     = mean(K2);
        JS3_mean    = mean(JS3);
        JK3_mean    = mean(JK3);

        if strcmp(mode,'time')
            if S1_mean < K2_mean
                wins_S1 = wins_S1 + 1;
            else
                wins_K2 = wins_K2 + 1;
            end

            if JS3_mean < JK3_mean
                wins_JS3 = wins_JS3 + 1;
            else
                wins_JK3 = wins_JK3 + 1;
            end
        else
            if S1_mean > K2_mean
                wins_S1 = wins_S1 + 1;
            else
                wins_K2 = wins_K2 + 1;
            end

            if JS3_mean > JK3_mean
                wins_JS3 = wins_JS3 + 1;
            else
                wins_JK3 = wins_JK3 + 1;
            end
        end
    end
    totalSolo   = wins_S1 + wins_K2;
    totalJoint  = wins_JS3 + wins_JK3;
    if wins_S1 > wins_K2
        winnerSolo = "S";
    elseif wins_K2 > wins_S1
        winnerSolo = "K";
    else
        winnerSolo = "tie";
    end

    if wins_JS3 > wins_JK3
        winnerJoint = "S";
    elseif wins_JK3 > wins_JS3
        winnerJoint = "K";
    else
        winnerJoint = "tie";
    end
    %% binomial test
    if totalSolo > 0
        kSolo = max(wins_S1, wins_K2);
        nSolo = totalSolo;

        p_oneSolo = binocdf(min(kSolo, nSolo-kSolo), nSolo, 0.5);

        p_binomSolo = min(1, 2 * p_oneSolo);
    else
        p_binomSolo = NaN;
    end

    signif_binomSolo = p_binomSolo <= 0.05;

    BinomialTableSolo.Direction(dirIdx) = dirIdx;
    BinomialTableSolo.Wins_S(dirIdx) = wins_S1;
    BinomialTableSolo.Wins_K(dirIdx) = wins_K2;
    BinomialTableSolo.Winner(dirIdx) = winnerSolo;
    BinomialTableSolo.p_binom(dirIdx) = p_binomSolo; 
    BinomialTableSolo.Significant(dirIdx) = signif_binomSolo;
    BinomialTableSolo.S_Percent(dirIdx) = wins_S1/totalSolo;
    BinomialTableSolo.K_Percent(dirIdx) = wins_K2/totalSolo;

    if totalJoint > 0
        kJoint = max(wins_JS3, wins_JK3);
        nJoint = totalJoint;

        p_oneJoint = binocdf(min(kJoint, nJoint-kJoint), nJoint, 0.5);

        p_binomJoint = min(1, 2 * p_oneJoint);
    else
        p_binomJoint = NaN;
    end

    signif_binomJoint = p_binomJoint <= 0.05;

    BinomialTableJoint.Direction(dirIdx) = dirIdx;
    BinomialTableJoint.Wins_S(dirIdx) = wins_JS3;
    BinomialTableJoint.Wins_K(dirIdx) = wins_JK3;
    BinomialTableJoint.Winner(dirIdx) = winnerJoint;
    BinomialTableJoint.p_binom(dirIdx) = p_binomJoint; 
    BinomialTableJoint.Significant(dirIdx) = signif_binomJoint;
    BinomialTableJoint.S_Percent(dirIdx) = wins_JS3/totalJoint;
    BinomialTableJoint.K_Percent(dirIdx) = wins_JK3/totalJoint;
    % winner (S o K)
    Data_result.Solo(1).Winner(1).(dirName) = winnerSolo;
    Data_result.Solo(1).S_wins(1).(dirName) = wins_S1;
    Data_result.Solo(1).K_wins(1).(dirName) = wins_K2;
    Data_result.Solo(1).total(1).(dirName)  = totalSolo;
    Data_result.Solo(1).S_perc(1).(dirName) = wins_S1/totalSolo;
    Data_result.Solo(1).K_perc(1).(dirName) = wins_K2/totalSolo;

    Data_result.Joint(1).Winner(1).(dirName) = winnerJoint;
    Data_result.Joint(1).S_wins(1).(dirName) = wins_JS3;
    Data_result.Joint(1).K_wins(1).(dirName) = wins_JK3;
    Data_result.Joint(1).total(1).(dirName)  = totalJoint;
    Data_result.Joint(1).S_perc(1).(dirName) = wins_JS3/totalJoint;
    Data_result.Joint(1).K_perc(1).(dirName) = wins_JK3/totalJoint;
end