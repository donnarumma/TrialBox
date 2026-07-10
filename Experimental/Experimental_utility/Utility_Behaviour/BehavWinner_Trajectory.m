function [Data_result, BinomialTableSolo, BinomialTableJoint, ComparisonTable] = BehavWinner_Trajectory(data_S, data_K, par)
num_dir = par.num_dir;
num_sessions = length(data_S);
inField = par.InField; % Parametro sotto test (es. 'exitTime')

% --- SETUP TABELLE BINOMIALI ---
VarNames = {'Metric', 'Direction','Wins_S','Wins_K','Winner','p_binom','Significant','S_Perc','K_Perc',...
    'Bias_S','Bias_K','Side_S','Side_K','SD_S','SD_K','Divergence','SocialShift_S','SocialShift_K'};
VarTypes = {'string', 'double','double','double','string','double','double','double','double',...
    'double','double','string','string','double','double','double','double','double'};

BinomialTableSolo = table('Size',[num_dir,length(VarNames)], 'VariableNames', VarNames, 'VariableTypes', VarTypes);
BinomialTableJoint = table('Size',[num_dir,length(VarNames)], 'VariableNames', VarNames, 'VariableTypes', VarTypes);

% --- SETUP COMPARISON TABLE ---
CompNames = {'Metric', 'Direction', 'Winner_Solo', 'Winner_Joint', 'BiasS_Solo', 'BiasS_Joint', 'Shift_S', 'Reliability_S',...
    'BiasK_Solo', 'BiasK_Joint', 'Shift_K', 'Reliability_K', 'Div_Solo', 'Div_Joint', 'Delta_Divergence'};
CompTypes = {'string', 'double', 'string', 'string', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', 'double'};
ComparisonTable = table('Size',[num_dir,length(CompNames)], 'VariableNames', CompNames, 'VariableTypes', CompTypes);

Data_result = struct();

for dirIdx = 1:num_dir
    dirName = sprintf('dir%d', dirIdx);
    wins_S1 = 0; wins_K2 = 0; wins_JS3 = 0; wins_JK3 = 0;

    % Preallocazione
    bS1_v = NaN(1, num_sessions); bK2_v = NaN(1, num_sessions);
    bS3_v = NaN(1, num_sessions); bK3_v = NaN(1, num_sessions);
    sdS1_v = NaN(1, num_sessions); sdK2_v = NaN(1, num_sessions);
    sdS3_v = NaN(1, num_sessions); sdK3_v = NaN(1, num_sessions);

    for nsess = 1:num_sessions
        % 1. Dati per il TEST (inField)
        vS1_raw = data_S(nsess).(inField)(1).(dirName);
        vK2_raw = data_K(nsess).(inField)(2).(dirName);
        vJS3_raw = data_S(nsess).(inField)(3).(dirName);
        vJK3_raw = data_K(nsess).(inField)(3).(dirName);

        % 2. Dati DESCRIZIONE (Sempre Angular Error per Bias/SD)
        angS1 = data_S(nsess).angular_error(1).(dirName);
        angK2 = data_K(nsess).angular_error(2).(dirName);
        angJS3 = data_S(nsess).angular_error(3).(dirName);
        angJK3 = data_K(nsess).angular_error(3).(dirName);

        bS1_v(nsess) = mean(angS1);   bK2_v(nsess) = mean(angK2);
        bS3_v(nsess) = mean(angJS3);  bK3_v(nsess) = mean(angJK3);
        sdS1_v(nsess) = std(angS1);   sdK2_v(nsess) = std(angK2);
        sdS3_v(nsess) = std(angJS3);  sdK3_v(nsess) = std(angJK3);

        valS1 = mean(vS1_raw); valK2 = mean(vK2_raw);
        valJS3 = mean(vJS3_raw); valJK3 = mean(vJK3_raw);


        % Assegnazione Wins
        if strcmp(inField, 'ExitCorr')
            if valS1 > valK2, wins_S1 = wins_S1 + 1; else, wins_K2 = wins_K2 + 1; end
            if valJS3 > valJK3, wins_JS3 = wins_JS3 + 1; else, wins_JK3 = wins_JK3 + 1; end
        else
            if valS1 < valK2, wins_S1 = wins_S1 + 1; else, wins_K2 = wins_K2 + 1; end
            if valJS3 < valJK3, wins_JS3 = wins_JS3 + 1; else, wins_JK3 = wins_JK3 + 1; end
        end
    end

    % Statistiche finali
    [wS, pS] = getBinom(wins_S1, wins_K2);
    [wJ, pJ] = getBinom(wins_JS3, wins_JK3);

    mBS1 = mean(bS1_v); mBK2 = mean(bK2_v);
    mBS3 = mean(bS3_v); mBK3 = mean(bK3_v);
    mSD1 = mean(sdS1_v); mSD2 = mean(sdK2_v); % SD baseline

    shiftS = abs(mBS3 - mBS1);
    shiftK = abs(mBK3 - mBK2);
    divSolo = abs(mBS1 - mBK2);
    divJoint = abs(mBS3 - mBK3);

    % --- POPOLAMENTO ---
    BinomialTableSolo(dirIdx,:) = {inField, dirIdx, wins_S1, wins_K2, wS, pS, pS<=0.05, ...
        100*(wins_S1/(wins_S1+wins_K2)), 100*(wins_K2/(wins_S1+wins_K2)), ...
        mBS1, mBK2, getSide(mBS1), getSide(mBK2), mSD1, mSD2, divSolo, 0, 0};

    BinomialTableJoint(dirIdx,:) = {inField, dirIdx, wins_JS3, wins_JK3, wJ, pJ, pJ<=0.05, ...
        100*(wins_JS3/(wins_JS3+wins_JK3)), 100*(wins_JK3/(wins_JS3+wins_JK3)), ...
        mBS3, mBK3, getSide(mBS3), getSide(mBK3), mean(sdS3_v), mean(sdK3_v), divJoint, shiftS, shiftK};

    % Reliability = Shift / SD_Solo (Valuta se lo spostamento supera il rumore motorio individuale)
    ComparisonTable(dirIdx,:) = {inField, dirIdx, wS, wJ, mBS1, mBS3, shiftS, shiftS/mSD1, ...
        mBK2, mBK3, shiftK, shiftK/mSD2, divSolo, divJoint, divJoint - divSolo};

    Data_result.Solo.Winner.(dirName) = wS;
    Data_result.Joint.Winner.(dirName) = wJ;
    Data_result.Comparison.DeltaDivergence.(dirName) = divJoint - divSolo;
end
end