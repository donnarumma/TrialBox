%% MAIN_SAPIENZA_BEHAVIOUR_NEURAL_RT_INCOMING_FOCUS
% Focused direction-wise RT analysis for incoming joint weights.
%
% The figure makes explicit whether S RT advantage is associated with
% S->K incoming-to-K weights (I_K) or K->S incoming-to-S weights (I_S).

clear;
close all;

analysisPms = struct();
analysisPms.inputRoot = '/TESTS/SAPIENZA/BEHAVIOUR/NeuralCorrelation';
analysisPms.outputRoot = fullfile(analysisPms.inputRoot, 'RTPattern');
analysisPms.visible = 'off';
analysisPms.directionIds = 1:8;
analysisPms.monkeySColor = [0.20 0.40 0.85];
analysisPms.monkeyKColor = [0.20 0.60 0.20];

plotDir = fullfile(analysisPms.outputRoot, 'plots');
pdfDir = fullfile(plotDir, 'PDF');
figDir = fullfile(plotDir, 'FIG');
if ~exist(pdfDir, 'dir')
    mkdir(pdfDir);
end
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

directionalPath = fullfile(analysisPms.inputRoot, ...
    'BehaviourNeural_directionalAllTargetCorrelations.csv');
globalPath = fullfile(analysisPms.inputRoot, ...
    'BehaviourNeural_correlations_sorted.csv');
if ~exist(directionalPath, 'file')
    error('%s:MissingInput', mfilename, 'Missing %s', directionalPath);
end
if ~exist(globalPath, 'file')
    error('%s:MissingInput', mfilename, 'Missing %s', globalPath);
end

directionalTable = readtable(directionalPath);
globalTable = readtable(globalPath);

focusPairs = cell2table({ ...
    'Frontal',  'M7', 'JointWeightSToK', 'I_K = S->K', 'Trajectory-aware Frontal M7'; ...
    'Frontal',  'M7', 'JointWeightKToS', 'I_S = K->S', 'Trajectory-aware Frontal M7'; ...
    'Parietal', 'M7', 'JointWeightSToK', 'I_K = S->K', 'Trajectory-aware Parietal M7'; ...
    'Parietal', 'M7', 'JointWeightKToS', 'I_S = K->S', 'Trajectory-aware Parietal M7'; ...
    'Parietal', 'M2', 'JointWeightSToK', 'I_K = S->K', 'Direct-coupling Parietal M2'; ...
    'Parietal', 'M2', 'JointWeightKToS', 'I_S = K->S', 'Direct-coupling Parietal M2'}, ...
    'VariableNames', {'Chamber', 'Model', 'NeuralReadout', 'ReadoutLabel', 'ModelLabel'});

focusSummary = local_build_focus_summary(directionalTable, globalTable, focusPairs);
writetable(focusSummary, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_RT_incomingWeightFocus.csv'));

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'RT advantage incoming-weight focus', ...
    'Position', [80 80 1250 1050]);
tileLayoutHandle = tiledlayout(3, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for focusPairCounter = 1:height(focusPairs)
    focusPair = focusPairs(focusPairCounter, :);
    nexttile;

    pairMask = strcmp(focusSummary.Chamber, focusPair.Chamber{1}) & ...
        strcmp(focusSummary.Model, focusPair.Model{1}) & ...
        strcmp(focusSummary.NeuralReadout, focusPair.NeuralReadout{1});
    pairTable = sortrows(focusSummary(pairMask, :), 'Direction');

    if strcmp(focusPair.NeuralReadout{1}, 'JointWeightSToK')
        barColor = analysisPms.monkeySColor;
    else
        barColor = analysisPms.monkeyKColor;
    end

    bar(pairTable.Direction, pairTable.DirectionRho, ...
        'FaceColor', barColor, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.78);
    hold on;
    yline(0, 'k-', 'LineWidth', 0.8);
    yline(pairTable.GlobalRho(1), '--', 'Color', [0.10 0.10 0.10], 'LineWidth', 1.3);
    grid on;
    xticks(analysisPms.directionIds);
    xticklabels(compose('D%d', analysisPms.directionIds));
    ylim([-0.75 0.75]);
    ylabel('Spearman rho');
    title({sprintf('%s %s', focusPair.Chamber{1}, focusPair.Model{1}), ...
        sprintf('%s; global rho=%+.2f, p=%s', focusPair.ReadoutLabel{1}, ...
        pairTable.GlobalRho(1), local_pvalue_text(pairTable.GlobalPValue(1)))}, ...
        'Interpreter', 'none');
end

title(tileLayoutHandle, ...
    'S RT advantage versus incoming joint weights by direction', ...
    'Interpreter', 'none');

figureStem = 'BehaviourNeural_RT_IncomingWeightFocus_DirectionalBars';
exportgraphics(figureHandle, fullfile(pdfDir, [figureStem '.pdf']), ...
    'ContentType', 'vector');
savefig(figureHandle, fullfile(figDir, [figureStem '.fig']));
close(figureHandle);

fprintf('Focused RT incoming-weight figure saved in %s\n', pdfDir);


function focusSummary = local_build_focus_summary(directionalTable, globalTable, focusPairs)
summaryRows = {};

for focusPairCounter = 1:height(focusPairs)
    focusPair = focusPairs(focusPairCounter, :);
    directionMask = strcmp(directionalTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(directionalTable.Model, focusPair.Model{1}) & ...
        strcmp(directionalTable.Predictor, 'SAdvantage_RT') & ...
        strcmp(directionalTable.Target, focusPair.NeuralReadout{1});
    pairDirectionalTable = sortrows(directionalTable(directionMask, :), 'Direction');

    globalMask = strcmp(globalTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(globalTable.Model, focusPair.Model{1}) & ...
        strcmp(globalTable.Predictor, 'SAdvantage_RT') & ...
        strcmp(globalTable.Target, focusPair.NeuralReadout{1});
    pairGlobalTable = globalTable(globalMask, :);

    if isempty(pairGlobalTable)
        globalRho = NaN;
        globalPValue = NaN;
        globalSampleCount = NaN;
    else
        globalRho = pairGlobalTable.Rho(1);
        globalPValue = pairGlobalTable.PValue(1);
        globalSampleCount = pairGlobalTable.N(1);
    end

    for directionCounter = 1:height(pairDirectionalTable)
        summaryRows(end + 1, :) = { ...
            focusPair.Chamber{1}, focusPair.Model{1}, focusPair.ModelLabel{1}, ...
            focusPair.NeuralReadout{1}, focusPair.ReadoutLabel{1}, ...
            pairDirectionalTable.Direction(directionCounter), ...
            pairDirectionalTable.Rho(directionCounter), ...
            pairDirectionalTable.PValue(directionCounter), ...
            pairDirectionalTable.N(directionCounter), ...
            globalRho, globalPValue, globalSampleCount}; %#ok<AGROW>
    end
end

focusSummary = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Model', 'ModelLabel', 'NeuralReadout', ...
    'ReadoutLabel', 'Direction', 'DirectionRho', 'DirectionPValue', ...
    'DirectionSampleCount', 'GlobalRho', 'GlobalPValue', 'GlobalSampleCount'});
end


function pValueText = local_pvalue_text(pValue)
if ~isfinite(pValue)
    pValueText = 'NA';
elseif pValue < 1e-3
    pValueText = sprintf('%.1e', pValue);
else
    pValueText = sprintf('%.3f', pValue);
end
end
