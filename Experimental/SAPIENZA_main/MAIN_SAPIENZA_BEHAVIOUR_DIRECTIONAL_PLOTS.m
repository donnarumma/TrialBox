%% MAIN_SAPIENZA_BEHAVIOUR_DIRECTIONAL_PLOTS
% Plot direction-resolved behavioural readouts from the RAW Sapienza output.
%
% The figure mirrors the neural PEB joint-action by-direction summaries:
%   left   = joint raw values for monkey S and monkey K;
%   middle = S advantage in joint action;
%   right  = joint benefit relative to each monkey's source-active baseline.

clear;
close all;

currentScriptDir = fileparts(mfilename('fullpath'));
trialBoxRoot = fileparts(fileparts(currentScriptDir));
addpath(fullfile(trialBoxRoot, 'Experimental', 'SAPIENZA_main'));

plotPms = struct();
plotPms.behaviourRoot = '/TESTS/SAPIENZA/BEHAVIOUR';
plotPms.outputRoot = fullfile(plotPms.behaviourRoot, 'DirectionalPlots');
plotPms.chamberNames = {'Frontal', 'Parietal'};
plotPms.fileConditionTag = 'D1D2D3D4D5D6D7D8_C1C2C3';
plotPms.metricNames = {'RT', 'PV', 'PVT', 'AE_abs', 'CMT', 'EC', 'ExT'};
plotPms.metricLabels = {'RT', 'PV', 'PVT', 'AE abs', 'CMT', 'EC', 'ExT'};
plotPms.directionIds = 1:8;
plotPms.directionLabels = {'D1 forward', 'D2 up-right', 'D3 right', ...
    'D4 down-right', 'D5 backward', 'D6 down-left', 'D7 left', 'D8 up-left'};
plotPms.directionAnglesDegrees = [0 45 90 135 180 225 270 315];
plotPms.monkeySColor = [0.20 0.40 0.85];
plotPms.monkeyKColor = [0.20 0.60 0.20];
plotPms.visible = 'off';

pdfDir = fullfile(plotPms.outputRoot, 'plots', 'PDF');
figDir = fullfile(plotPms.outputRoot, 'plots', 'FIG');
if ~exist(pdfDir, 'dir')
    mkdir(pdfDir);
end
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

summaryRows = {};

for chamberCounter = 1:numel(plotPms.chamberNames)
    chamberName = plotPms.chamberNames{chamberCounter};
    chamberPath = fullfile(plotPms.behaviourRoot, chamberName, ...
        sprintf('Behaviour_RAW_%s_AllSessions_%s.mat', chamberName, plotPms.fileConditionTag));

    if ~exist(chamberPath, 'file')
        warning('%s missing behaviour file: %s', mfilename, chamberPath);
        continue;
    end

    loadedBehaviour = load(chamberPath, 'Data');
    chamberBehaviour = loadedBehaviour.Data.Behav;
    sessionNames = fieldnames(chamberBehaviour);

    fprintf('Plotting directional behaviour for %s (%d sessions).\n', ...
        chamberName, numel(sessionNames));

    [figureHandle, chamberSummaryRows] = local_plot_chamber_directional_behaviour( ...
        chamberBehaviour, sessionNames, chamberName, plotPms);

    summaryRows = [summaryRows; chamberSummaryRows]; %#ok<AGROW>

    figureStem = sprintf('BehaviourDirectional_%s_byDirection', chamberName);
    exportgraphics(figureHandle, fullfile(pdfDir, [figureStem '.pdf']), ...
        'ContentType', 'vector');
    savefig(figureHandle, fullfile(figDir, [figureStem '.fig']));
    close(figureHandle);
end

summaryTable = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Metric', 'Direction', 'N', ...
    'MeanSAct', 'SemSAct', 'MeanKAct', 'SemKAct', ...
    'MeanSJoint', 'SemSJoint', 'MeanKJoint', 'SemKJoint', ...
    'MeanSAdvantage', 'SemSAdvantage', ...
    'MeanSJointBenefit', 'SemSJointBenefit', ...
    'MeanKJointBenefit', 'SemKJointBenefit'});

writetable(summaryTable, fullfile(plotPms.outputRoot, 'BehaviourDirectional_summary.csv'));
save(fullfile(plotPms.outputRoot, 'BehaviourDirectional_summary.mat'), ...
    'plotPms', 'summaryTable');

[preferenceFigureHandle, preferenceTable] = local_plot_directional_preference(summaryTable, plotPms);
preferenceStem = 'BehaviourDirectional_LeaderFollowerPreference';
writetable(preferenceTable, fullfile(plotPms.outputRoot, ...
    [preferenceStem '.csv']));
exportgraphics(preferenceFigureHandle, fullfile(pdfDir, [preferenceStem '.pdf']), ...
    'ContentType', 'vector');
savefig(preferenceFigureHandle, fullfile(figDir, [preferenceStem '.fig']));
close(preferenceFigureHandle);

[aeFigureHandle, aeTable] = local_plot_ae_abs_condition_comparison(summaryTable, plotPms);
aeFigureStem = 'BehaviourDirectional_AEAbs_SourceActiveVsJoint';
writetable(aeTable, fullfile(plotPms.outputRoot, [aeFigureStem '.csv']));
exportgraphics(aeFigureHandle, fullfile(pdfDir, [aeFigureStem '.pdf']), ...
    'ContentType', 'vector');
savefig(aeFigureHandle, fullfile(figDir, [aeFigureStem '.fig']));
close(aeFigureHandle);

[aeKeyFigureHandle, aeKeyTable] = local_plot_ae_abs_key_pattern(summaryTable, plotPms);
aeKeyFigureStem = 'BehaviourDirectional_AEAbs_KeyPattern';
writetable(aeKeyTable, fullfile(plotPms.outputRoot, [aeKeyFigureStem '.csv']));
exportgraphics(aeKeyFigureHandle, fullfile(pdfDir, [aeKeyFigureStem '.pdf']), ...
    'ContentType', 'vector');
savefig(aeKeyFigureHandle, fullfile(figDir, [aeKeyFigureStem '.fig']));
close(aeKeyFigureHandle);

[aePolarFigureHandle, aePolarTable] = local_plot_ae_abs_polar_pattern(aeKeyTable, plotPms);
aePolarFigureStem = 'BehaviourDirectional_AEAbs_PolarPattern';
writetable(aePolarTable, fullfile(plotPms.outputRoot, [aePolarFigureStem '.csv']));
exportgraphics(aePolarFigureHandle, fullfile(pdfDir, [aePolarFigureStem '.pdf']), ...
    'ContentType', 'vector');
savefig(aePolarFigureHandle, fullfile(figDir, [aePolarFigureStem '.fig']));
close(aePolarFigureHandle);

fprintf('Directional behaviour plots saved in %s\n', plotPms.outputRoot);


function [figureHandle, summaryRows] = local_plot_chamber_directional_behaviour( ...
    chamberBehaviour, sessionNames, chamberName, plotPms)

metricCount = numel(plotPms.metricNames);
directionIds = plotPms.directionIds;
directionLabels = compose('D%d', directionIds);
summaryRows = {};

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', plotPms.visible, ...
    'Name', sprintf('%s raw behaviour by direction', chamberName), ...
    'Position', [80 80 1480 max(850, 210 * metricCount)]);

tileLayoutHandle = tiledlayout(metricCount, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');
title(tileLayoutHandle, ...
    sprintf('%s: raw behavioural features by direction', chamberName), ...
    'Interpreter', 'none');

for metricCounter = 1:metricCount
    metricName = plotPms.metricNames{metricCounter};
    metricLabel = plotPms.metricLabels{metricCounter};
    metricData = local_collect_metric_by_direction( ...
        chamberBehaviour, sessionNames, metricName, directionIds);

    [meanSJoint, semSJoint, sampleCount] = local_mean_sem(metricData.sJoint);
    [meanKJoint, semKJoint] = local_mean_sem(metricData.kJoint);
    [meanSAct, semSAct] = local_mean_sem(metricData.sAct);
    [meanKAct, semKAct] = local_mean_sem(metricData.kAct);
    [meanSAdvantage, semSAdvantage] = local_mean_sem(metricData.sAdvantage);
    [meanSJointBenefit, semSJointBenefit] = local_mean_sem(metricData.sJointBenefit);
    [meanKJointBenefit, semKJointBenefit] = local_mean_sem(metricData.kJointBenefit);

    for directionCounter = 1:numel(directionIds)
        summaryRows(end + 1, :) = { ...
            chamberName, metricName, directionIds(directionCounter), ...
            sampleCount(directionCounter), ...
            meanSAct(directionCounter), semSAct(directionCounter), ...
            meanKAct(directionCounter), semKAct(directionCounter), ...
            meanSJoint(directionCounter), semSJoint(directionCounter), ...
            meanKJoint(directionCounter), semKJoint(directionCounter), ...
            meanSAdvantage(directionCounter), semSAdvantage(directionCounter), ...
            meanSJointBenefit(directionCounter), semSJointBenefit(directionCounter), ...
            meanKJointBenefit(directionCounter), semKJointBenefit(directionCounter)}; %#ok<AGROW>
    end

    nexttile;
    local_plot_two_series(directionIds, meanSJoint, semSJoint, meanKJoint, semKJoint, ...
        plotPms.monkeySColor, plotPms.monkeyKColor);
    if metricCounter == 1
        title('Joint raw values');
        legend({'S joint', 'K joint'}, 'Location', 'best');
    end
    ylabel(metricLabel, 'Interpreter', 'none');
    local_format_direction_axis(directionIds, directionLabels);

    nexttile;
    local_plot_signed_advantage(directionIds, meanSAdvantage, semSAdvantage, ...
        plotPms.monkeySColor, plotPms.monkeyKColor);
    yline(0, 'k-');
    if metricCounter == 1
        title('S advantage in joint');
        legend({'S advantage', 'K advantage'}, 'Location', 'best');
    end
    ylabel(sprintf('%s advantage', metricLabel), 'Interpreter', 'none');
    local_format_direction_axis(directionIds, directionLabels);

    nexttile;
    local_plot_two_series(directionIds, meanSJointBenefit, semSJointBenefit, ...
        meanKJointBenefit, semKJointBenefit, plotPms.monkeySColor, plotPms.monkeyKColor);
    yline(0, 'k-');
    if metricCounter == 1
        title('Joint benefit vs source-active baseline');
        legend({'S benefit', 'K benefit'}, 'Location', 'best');
    end
    ylabel(sprintf('%s benefit', metricLabel), 'Interpreter', 'none');
    local_format_direction_axis(directionIds, directionLabels);
end
end


function metricData = local_collect_metric_by_direction( ...
    chamberBehaviour, sessionNames, metricName, directionIds)

sessionCount = numel(sessionNames);
directionCount = numel(directionIds);
metricData = struct();
metricData.sAct = nan(sessionCount, directionCount);
metricData.kAct = nan(sessionCount, directionCount);
metricData.sJoint = nan(sessionCount, directionCount);
metricData.kJoint = nan(sessionCount, directionCount);
metricData.sAdvantage = nan(sessionCount, directionCount);
metricData.sJointBenefit = nan(sessionCount, directionCount);
metricData.kJointBenefit = nan(sessionCount, directionCount);

for sessionCounter = 1:sessionCount
    sessionName = sessionNames{sessionCounter};
    sessionBehaviour = chamberBehaviour.(sessionName);

    for directionCounter = 1:directionCount
        directionId = directionIds(directionCounter);
        sActIndex = directionId;
        kActIndex = 8 + directionId;
        jointIndex = 16 + directionId;

        metricData.sAct(sessionCounter, directionCounter) = ...
            local_get_metric_value(sessionBehaviour, 'S', sActIndex, metricName);
        metricData.kAct(sessionCounter, directionCounter) = ...
            local_get_metric_value(sessionBehaviour, 'K', kActIndex, metricName);
        metricData.sJoint(sessionCounter, directionCounter) = ...
            local_get_metric_value(sessionBehaviour, 'S', jointIndex, metricName);
        metricData.kJoint(sessionCounter, directionCounter) = ...
            local_get_metric_value(sessionBehaviour, 'K', jointIndex, metricName);

        metricData.sAdvantage(sessionCounter, directionCounter) = ...
            local_s_advantage(metricName, ...
            metricData.sJoint(sessionCounter, directionCounter), ...
            metricData.kJoint(sessionCounter, directionCounter));
        metricData.sJointBenefit(sessionCounter, directionCounter) = ...
            local_joint_benefit(metricName, ...
            metricData.sAct(sessionCounter, directionCounter), ...
            metricData.sJoint(sessionCounter, directionCounter));
        metricData.kJointBenefit(sessionCounter, directionCounter) = ...
            local_joint_benefit(metricName, ...
            metricData.kAct(sessionCounter, directionCounter), ...
            metricData.kJoint(sessionCounter, directionCounter));
    end
end
end


function metricValue = local_get_metric_value(sessionBehaviour, monkeyName, behaviourIndex, metricName)
metricValue = NaN;
if ~isfield(sessionBehaviour, monkeyName)
    return;
end

monkeyBehaviour = sessionBehaviour.(monkeyName);
if numel(monkeyBehaviour) < behaviourIndex
    return;
end
if ~isfield(monkeyBehaviour(behaviourIndex), metricName)
    return;
end

metricValue = monkeyBehaviour(behaviourIndex).(metricName);
end


function sAdvantage = local_s_advantage(metricName, sJointValue, kJointValue)
if local_is_lower_better(metricName)
    sAdvantage = kJointValue - sJointValue;
else
    sAdvantage = sJointValue - kJointValue;
end
end


function jointBenefit = local_joint_benefit(metricName, actValue, jointValue)
if local_is_lower_better(metricName)
    jointBenefit = actValue - jointValue;
else
    jointBenefit = jointValue - actValue;
end
end


function isLowerBetter = local_is_lower_better(metricName)
isLowerBetter = ismember(metricName, {'RT', 'PVT', 'AE_abs', 'CMT', 'ExT'});
end


function [meanValues, semValues, sampleCounts] = local_mean_sem(valueMatrix)
sampleCounts = sum(isfinite(valueMatrix), 1);
meanValues = mean(valueMatrix, 1, 'omitnan');
standardDeviation = std(valueMatrix, 0, 1, 'omitnan');
semValues = standardDeviation ./ sqrt(sampleCounts);
semValues(sampleCounts == 0) = NaN;
end


function local_plot_two_series(directionIds, meanFirst, semFirst, meanSecond, semSecond, ...
    firstColor, secondColor)
hold on;
errorbar(directionIds, meanFirst, semFirst, '-o', ...
    'Color', firstColor, ...
    'MarkerFaceColor', firstColor, ...
    'LineWidth', 1.5, ...
    'CapSize', 4);
errorbar(directionIds, meanSecond, semSecond, '-o', ...
    'Color', secondColor, ...
    'MarkerFaceColor', secondColor, ...
    'LineWidth', 1.5, ...
    'CapSize', 4);
grid on;
end


function local_plot_signed_advantage(directionIds, meanValues, semValues, positiveColor, negativeColor)
hold on;

plot(nan, nan, '-o', ...
    'Color', positiveColor, ...
    'MarkerFaceColor', positiveColor, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'S advantage');
plot(nan, nan, '-o', ...
    'Color', negativeColor, ...
    'MarkerFaceColor', negativeColor, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'K advantage');

for segmentCounter = 1:(numel(directionIds) - 1)
    segmentDirections = directionIds(segmentCounter:(segmentCounter + 1));
    segmentValues = meanValues(segmentCounter:(segmentCounter + 1));
    if any(~isfinite(segmentValues))
        continue;
    end

    if mean(segmentValues, 'omitnan') >= 0
        segmentColor = positiveColor;
    else
        segmentColor = negativeColor;
    end

    plot(segmentDirections, segmentValues, '-', ...
        'Color', segmentColor, ...
        'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
end

for pointCounter = 1:numel(directionIds)
    if ~isfinite(meanValues(pointCounter))
        continue;
    end

    if meanValues(pointCounter) >= 0
        pointColor = positiveColor;
    else
        pointColor = negativeColor;
    end

    errorbar(directionIds(pointCounter), meanValues(pointCounter), semValues(pointCounter), ...
        'o', ...
        'Color', pointColor, ...
        'MarkerFaceColor', pointColor, ...
        'MarkerEdgeColor', pointColor, ...
        'LineStyle', 'none', ...
        'LineWidth', 1.5, ...
        'CapSize', 4, ...
        'HandleVisibility', 'off');
end

grid on;
end


function [figureHandle, preferenceTable] = local_plot_directional_preference(summaryTable, plotPms)
preferenceRows = {};
directionIds = plotPms.directionIds;
directionLabels = compose('D%d', directionIds);
metricNames = plotPms.metricNames;
metricLabels = plotPms.metricLabels;

for rowCounter = 1:height(summaryTable)
    jointAdvantage = summaryTable.MeanSAdvantage(rowCounter);
    benefitBalance = summaryTable.MeanSJointBenefit(rowCounter) - ...
        summaryTable.MeanKJointBenefit(rowCounter);

    if jointAdvantage >= 0
        jointLeader = 'S';
    else
        jointLeader = 'K';
    end

    if benefitBalance >= 0
        benefitLeader = 'S';
    else
        benefitLeader = 'K';
    end

    preferenceRows(end + 1, :) = { ...
        summaryTable.Chamber{rowCounter}, summaryTable.Metric{rowCounter}, ...
        summaryTable.Direction(rowCounter), jointAdvantage, jointLeader, ...
        benefitBalance, benefitLeader}; %#ok<AGROW>
end

preferenceTable = cell2table(preferenceRows, ...
    'VariableNames', {'Chamber', 'Metric', 'Direction', ...
    'JointAdvantage', 'JointLeader', 'JointBenefitBalance', 'BenefitLeader'});

maxAbsValue = max(abs([preferenceTable.JointAdvantage; preferenceTable.JointBenefitBalance]), [], 'omitnan');
if ~isfinite(maxAbsValue) || maxAbsValue == 0
    maxAbsValue = 1;
end

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', plotPms.visible, ...
    'Name', 'Directional behavioural leader-follower preference', ...
    'Position', [80 80 1320 840]);

tiledlayout(numel(plotPms.chamberNames), 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for chamberCounter = 1:numel(plotPms.chamberNames)
    chamberName = plotPms.chamberNames{chamberCounter};
    chamberMask = strcmp(preferenceTable.Chamber, chamberName);

    jointAdvantageMatrix = local_preference_matrix( ...
        preferenceTable(chamberMask, :), metricNames, directionIds, 'JointAdvantage');
    benefitBalanceMatrix = local_preference_matrix( ...
        preferenceTable(chamberMask, :), metricNames, directionIds, 'JointBenefitBalance');

    nexttile;
    local_plot_preference_heatmap(jointAdvantageMatrix, metricLabels, directionLabels, ...
        maxAbsValue, plotPms, sprintf('%s: joint performance leader', chamberName), ...
        'A_S^{(m)}');

    nexttile;
    local_plot_preference_heatmap(benefitBalanceMatrix, metricLabels, directionLabels, ...
        maxAbsValue, plotPms, sprintf('%s: joint-benefit leader', chamberName), ...
        'B_S^{(m)} - B_K^{(m)}');
end

sgtitle('Directional behavioural preference: S-positive, K-negative', ...
    'Interpreter', 'tex');
end


function valueMatrix = local_preference_matrix(preferenceTable, metricNames, directionIds, valueName)
valueMatrix = nan(numel(metricNames), numel(directionIds));

for metricCounter = 1:numel(metricNames)
    metricName = metricNames{metricCounter};
    for directionCounter = 1:numel(directionIds)
        directionId = directionIds(directionCounter);
        rowMask = strcmp(preferenceTable.Metric, metricName) & ...
            preferenceTable.Direction == directionId;
        if any(rowMask)
            valueMatrix(metricCounter, directionCounter) = ...
                preferenceTable.(valueName)(find(rowMask, 1, 'first'));
        end
    end
end
end


function local_plot_preference_heatmap(valueMatrix, metricLabels, directionLabels, ...
    maxAbsValue, plotPms, titleText, colorbarLabel)
imageHandle = imagesc(valueMatrix, [-maxAbsValue maxAbsValue]);
set(imageHandle, 'AlphaData', isfinite(valueMatrix));
set(gca, 'Color', [0.93 0.93 0.93]);
colormap(gca, local_s_k_colormap(256, plotPms.monkeyKColor, plotPms.monkeySColor));
colorbarHandle = colorbar;
ylabel(colorbarHandle, colorbarLabel, 'Interpreter', 'tex');
xticks(1:numel(directionLabels));
xticklabels(directionLabels);
yticks(1:numel(metricLabels));
yticklabels(metricLabels);
xlabel('Direction');
title(titleText, 'Interpreter', 'tex');
grid on;

for metricCounter = 1:size(valueMatrix, 1)
    for directionCounter = 1:size(valueMatrix, 2)
        preferenceValue = valueMatrix(metricCounter, directionCounter);
        if ~isfinite(preferenceValue)
            continue;
        end

        if preferenceValue >= 0
            preferenceLabel = 'S';
        else
            preferenceLabel = 'K';
        end

        if abs(preferenceValue) > 0.55 * maxAbsValue
            textColor = [1 1 1];
        else
            textColor = [0.05 0.05 0.05];
        end

        text(directionCounter, metricCounter, preferenceLabel, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontWeight', 'bold', ...
            'Color', textColor);
    end
end
end


function colorMap = local_s_k_colormap(colorCount, monkeyKColor, monkeySColor)
halfColorCount = colorCount / 2;
kToWhite = [linspace(monkeyKColor(1), 1, halfColorCount)', ...
    linspace(monkeyKColor(2), 1, halfColorCount)', ...
    linspace(monkeyKColor(3), 1, halfColorCount)'];
whiteToS = [linspace(1, monkeySColor(1), halfColorCount)', ...
    linspace(1, monkeySColor(2), halfColorCount)', ...
    linspace(1, monkeySColor(3), halfColorCount)'];
colorMap = [kToWhite; whiteToS];
end


function [figureHandle, aeTable] = local_plot_ae_abs_condition_comparison(summaryTable, plotPms)
aeTable = summaryTable(strcmp(summaryTable.Metric, 'AE_abs'), :);
if isempty(aeTable)
    figureHandle = figure('Visible', plotPms.visible);
    return;
end

directionIds = plotPms.directionIds;
directionLabels = compose('D%d', directionIds);
lightSColor = local_lighten_color(plotPms.monkeySColor, 0.58);
lightKColor = local_lighten_color(plotPms.monkeyKColor, 0.58);

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', plotPms.visible, ...
    'Name', 'AE abs source-active versus joint by direction', ...
    'Position', [80 80 1460 760]);
tiledlayout(numel(plotPms.chamberNames), 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for chamberCounter = 1:numel(plotPms.chamberNames)
    chamberName = plotPms.chamberNames{chamberCounter};
    chamberTable = aeTable(strcmp(aeTable.Chamber, chamberName), :);
    chamberTable = sortrows(chamberTable, 'Direction');

    nexttile;
    local_plot_condition_pair(directionIds, ...
        chamberTable.MeanSAct, chamberTable.SemSAct, ...
        chamberTable.MeanSJoint, chamberTable.SemSJoint, ...
        lightSColor, plotPms.monkeySColor, ...
        {'S ActS,ObsK', 'S Joint'});
    title(sprintf('%s S: source-active vs joint', chamberName), ...
        'Interpreter', 'none');
    ylabel('|AE| (lower is better)');
    local_format_direction_axis(directionIds, directionLabels);

    nexttile;
    local_plot_condition_pair(directionIds, ...
        chamberTable.MeanKAct, chamberTable.SemKAct, ...
        chamberTable.MeanKJoint, chamberTable.SemKJoint, ...
        lightKColor, plotPms.monkeyKColor, ...
        {'K ObsS,ActK', 'K Joint'});
    title(sprintf('%s K: source-active vs joint', chamberName), ...
        'Interpreter', 'none');
    ylabel('|AE| (lower is better)');
    local_format_direction_axis(directionIds, directionLabels);

    nexttile;
    local_plot_two_series(directionIds, ...
        chamberTable.MeanSJointBenefit, chamberTable.SemSJointBenefit, ...
        chamberTable.MeanKJointBenefit, chamberTable.SemKJointBenefit, ...
        plotPms.monkeySColor, plotPms.monkeyKColor);
    yline(0, 'k-');
    title(sprintf('%s: joint benefit in |AE|', chamberName), ...
        'Interpreter', 'none');
    ylabel('Benefit = source-active - joint');
    if chamberCounter == 1
        legend({'S benefit', 'K benefit'}, 'Location', 'best');
    end
    local_format_direction_axis(directionIds, directionLabels);
end

sgtitle('Absolute angular error by condition and direction', ...
    'Interpreter', 'none');
end


function [figureHandle, aeKeyTable] = local_plot_ae_abs_key_pattern(summaryTable, plotPms)
aeTable = summaryTable(strcmp(summaryTable.Metric, 'AE_abs'), :);
if isempty(aeTable)
    figureHandle = figure('Visible', plotPms.visible);
    aeKeyTable = table();
    return;
end

directionIds = plotPms.directionIds(:);
directionLabels = plotPms.directionLabels(:);
shortDirectionLabels = compose('D%d', directionIds);
lightSColor = local_lighten_color(plotPms.monkeySColor, 0.58);
lightKColor = local_lighten_color(plotPms.monkeyKColor, 0.58);

[meanSAct, semSAct, pooledCount] = local_pool_direction_stats( ...
    aeTable, directionIds, 'MeanSAct', 'SemSAct');
[meanSJoint, semSJoint] = local_pool_direction_stats( ...
    aeTable, directionIds, 'MeanSJoint', 'SemSJoint');
[meanKAct, semKAct] = local_pool_direction_stats( ...
    aeTable, directionIds, 'MeanKAct', 'SemKAct');
[meanKJoint, semKJoint] = local_pool_direction_stats( ...
    aeTable, directionIds, 'MeanKJoint', 'SemKJoint');
[meanSJointBenefit, semSJointBenefit] = local_pool_direction_stats( ...
    aeTable, directionIds, 'MeanSJointBenefit', 'SemSJointBenefit');
[meanKJointBenefit, semKJointBenefit] = local_pool_direction_stats( ...
    aeTable, directionIds, 'MeanKJointBenefit', 'SemKJointBenefit');

aeKeyTable = table(directionIds, directionLabels, pooledCount, ...
    meanSAct, semSAct, meanSJoint, semSJoint, ...
    meanKAct, semKAct, meanKJoint, semKJoint, ...
    meanSJointBenefit, semSJointBenefit, ...
    meanKJointBenefit, semKJointBenefit, ...
    'VariableNames', {'Direction', 'DirectionLabel', 'N', ...
    'MeanSAct', 'SemSAct', 'MeanSJoint', 'SemSJoint', ...
    'MeanKAct', 'SemKAct', 'MeanKJoint', 'SemKJoint', ...
    'MeanSJointBenefit', 'SemSJointBenefit', ...
    'MeanKJointBenefit', 'SemKJointBenefit'});

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', plotPms.visible, ...
    'Name', 'AE abs key behavioural pattern', ...
    'Position', [80 80 1440 900]);
tileLayoutHandle = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;
local_plot_condition_pair(directionIds, meanSAct, semSAct, ...
    meanSJoint, semSJoint, lightSColor, plotPms.monkeySColor, ...
    {'S ActS,ObsK', 'S Joint'});
local_highlight_direction_band([1 2], [0.93 0.96 1.00]);
local_format_direction_axis(directionIds, shortDirectionLabels);
title('S: forward/up-right errors are selectively reduced in joint', ...
    'Interpreter', 'none');
ylabel('|AE| (deg, lower is better)');
local_annotate_direction_group(2.05, max([meanSAct; meanSJoint], [], 'omitnan'), ...
    {'D1-D2: high S error', 'joint correction'});

nexttile;
local_plot_condition_pair(directionIds, meanKAct, semKAct, ...
    meanKJoint, semKJoint, lightKColor, plotPms.monkeyKColor, ...
    {'K ObsS,ActK', 'K Joint'});
local_highlight_direction_band([8 8], [0.94 1.00 0.94]);
local_format_direction_axis(directionIds, shortDirectionLabels);
title('K: up-left remains problematic in joint', ...
    'Interpreter', 'none');
ylabel('|AE| (deg, lower is better)');
local_annotate_direction_group(8, max([meanKAct; meanKJoint], [], 'omitnan'), ...
    {'D8: no', 'joint rescue'});

nexttile;
local_plot_benefit_bars(directionIds, meanSJointBenefit, semSJointBenefit, ...
    plotPms.monkeySColor);
local_highlight_direction_band([1 2], [0.93 0.96 1.00]);
local_format_direction_axis(directionIds, shortDirectionLabels);
title('S joint benefit: positive means lower |AE| in joint', ...
    'Interpreter', 'none');
ylabel('S source-active - S joint');

nexttile;
local_plot_benefit_bars(directionIds, meanKJointBenefit, semKJointBenefit, ...
    plotPms.monkeyKColor);
local_highlight_direction_band([8 8], [0.94 1.00 0.94]);
local_format_direction_axis(directionIds, shortDirectionLabels);
title('K joint benefit: D8 worsens in joint', ...
    'Interpreter', 'none');
ylabel('K source-active - K joint');

title(tileLayoutHandle, ...
    'Absolute angular error: direction-specific joint correction and failure modes', ...
    'Interpreter', 'none');
end


function [figureHandle, aePolarTable] = local_plot_ae_abs_polar_pattern(aeKeyTable, plotPms)
if isempty(aeKeyTable)
    figureHandle = figure('Visible', plotPms.visible);
    aePolarTable = table();
    return;
end

aePolarTable = aeKeyTable;
thetaRadians = deg2rad([plotPms.directionAnglesDegrees(:); 360]);
thetaTickLabels = compose('D%d', plotPms.directionIds);
lightSColor = local_lighten_color(plotPms.monkeySColor, 0.58);
lightKColor = local_lighten_color(plotPms.monkeyKColor, 0.58);
maxRadius = max([aeKeyTable.MeanSAct; aeKeyTable.MeanSJoint; ...
    aeKeyTable.MeanKAct; aeKeyTable.MeanKJoint], [], 'omitnan');
maxRadius = ceil(maxRadius / 10) * 10;

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', plotPms.visible, ...
    'Name', 'AE abs polar behavioural pattern', ...
    'Position', [80 80 1350 720]);
tileLayoutHandle = tiledlayout(1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

sPolarAxes = polaraxes(tileLayoutHandle);
sPolarAxes.Layout.Tile = 1;
local_plot_ae_abs_polar_axes(sPolarAxes, thetaRadians, ...
    aeKeyTable.MeanSAct, aeKeyTable.MeanSJoint, ...
    lightSColor, plotPms.monkeySColor, maxRadius, thetaTickLabels, ...
    'S: |AE| compass', {'S ActS,ObsK', 'S Joint'});

kPolarAxes = polaraxes(tileLayoutHandle);
kPolarAxes.Layout.Tile = 2;
local_plot_ae_abs_polar_axes(kPolarAxes, thetaRadians, ...
    aeKeyTable.MeanKAct, aeKeyTable.MeanKJoint, ...
    lightKColor, plotPms.monkeyKColor, maxRadius, thetaTickLabels, ...
    'K: |AE| compass', {'K ObsS,ActK', 'K Joint'});

title(tileLayoutHandle, ...
    'Polar view of absolute angular error by movement direction', ...
    'Interpreter', 'none');
end


function local_plot_ae_abs_polar_axes(polarAxesHandle, thetaRadians, ...
    sourceActiveValues, jointValues, sourceActiveColor, jointColor, ...
    maxRadius, thetaTickLabels, plotTitle, legendLabels)
sourceActiveValues = [sourceActiveValues(:); sourceActiveValues(1)];
jointValues = [jointValues(:); jointValues(1)];

hold(polarAxesHandle, 'on');
polarplot(polarAxesHandle, thetaRadians, sourceActiveValues, '--o', ...
    'Color', sourceActiveColor, ...
    'MarkerFaceColor', sourceActiveColor, ...
    'LineWidth', 1.6, ...
    'MarkerSize', 6);
polarplot(polarAxesHandle, thetaRadians, jointValues, '-o', ...
    'Color', jointColor, ...
    'MarkerFaceColor', jointColor, ...
    'LineWidth', 2.2, ...
    'MarkerSize', 6);
polarAxesHandle.ThetaZeroLocation = 'top';
polarAxesHandle.ThetaDir = 'clockwise';
polarAxesHandle.ThetaTick = 0:45:315;
polarAxesHandle.ThetaTickLabel = thetaTickLabels;
polarAxesHandle.RLim = [0 maxRadius];
polarAxesHandle.RAxis.Label.String = '|AE| (deg)';
polarAxesHandle.FontSize = 10;
title(polarAxesHandle, plotTitle, ...
    'Interpreter', 'none', ...
    'FontWeight', 'bold');
legend(polarAxesHandle, legendLabels, ...
    'Location', 'southoutside', ...
    'Interpreter', 'none');
end


function [pooledMeans, pooledSems, pooledCounts] = local_pool_direction_stats( ...
    summaryTable, directionIds, meanColumnName, semColumnName)
pooledMeans = nan(numel(directionIds), 1);
pooledSems = nan(numel(directionIds), 1);
pooledCounts = zeros(numel(directionIds), 1);

for directionCounter = 1:numel(directionIds)
    directionMask = summaryTable.Direction == directionIds(directionCounter);
    directionRows = summaryTable(directionMask, :);
    groupCounts = double(directionRows.N);
    groupMeans = directionRows.(meanColumnName);
    groupSems = directionRows.(semColumnName);
    validMask = groupCounts > 0 & ~isnan(groupMeans) & ~isnan(groupSems);
    groupCounts = groupCounts(validMask);
    groupMeans = groupMeans(validMask);
    groupSems = groupSems(validMask);

    if isempty(groupCounts)
        continue;
    end

    totalCount = sum(groupCounts);
    pooledCounts(directionCounter) = totalCount;
    pooledMeans(directionCounter) = sum(groupCounts .* groupMeans) ./ totalCount;

    if totalCount > 1
        groupVariances = (groupSems .^ 2) .* groupCounts;
        withinSquares = sum((groupCounts - 1) .* groupVariances);
        betweenSquares = sum(groupCounts .* ...
            (groupMeans - pooledMeans(directionCounter)) .^ 2);
        pooledVariance = (withinSquares + betweenSquares) ./ (totalCount - 1);
        pooledSems(directionCounter) = sqrt(pooledVariance ./ totalCount);
    end
end
end


function local_plot_benefit_bars(directionIds, meanBenefit, semBenefit, barColor)
hold on;
bar(directionIds, meanBenefit, ...
    'FaceColor', barColor, ...
    'FaceAlpha', 0.78, ...
    'EdgeColor', 'none');
errorbar(directionIds, meanBenefit, semBenefit, 'k.', ...
    'LineWidth', 1.2, ...
    'CapSize', 4);
yline(0, 'k-');
grid on;
end


function local_highlight_direction_band(directionRange, bandColor)
axisHandle = gca;
axisLimits = axisHandle.YLim;
patch([directionRange(1) - 0.5, directionRange(end) + 0.5, ...
    directionRange(end) + 0.5, directionRange(1) - 0.5], ...
    [axisLimits(1), axisLimits(1), axisLimits(2), axisLimits(2)], ...
    bandColor, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.32, ...
    'HandleVisibility', 'off');
uistack(findobj(axisHandle, 'Type', 'patch'), 'bottom');
end


function local_annotate_direction_group(xCoordinate, yReference, annotationText)
axisHandle = gca;
axisLimits = axisHandle.YLim;
yCoordinate = axisLimits(1) + 0.88 .* (axisLimits(2) - axisLimits(1));
if ~isnan(yReference)
    yCoordinate = min(yCoordinate, yReference .* 0.96);
end
text(xCoordinate, yCoordinate, annotationText, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 9, ...
    'FontWeight', 'bold', ...
    'BackgroundColor', 'w', ...
    'Margin', 2, ...
    'Interpreter', 'none');
end


function local_plot_condition_pair(directionIds, meanSourceActive, semSourceActive, ...
    meanJoint, semJoint, sourceActiveColor, jointColor, legendLabels)
hold on;
errorbar(directionIds, meanSourceActive, semSourceActive, '--o', ...
    'Color', sourceActiveColor, ...
    'MarkerFaceColor', sourceActiveColor, ...
    'LineWidth', 1.4, ...
    'CapSize', 4);
errorbar(directionIds, meanJoint, semJoint, '-o', ...
    'Color', jointColor, ...
    'MarkerFaceColor', jointColor, ...
    'LineWidth', 1.8, ...
    'CapSize', 4);
grid on;
legend(legendLabels, 'Location', 'best');
end


function lightColor = local_lighten_color(baseColor, blendFraction)
lightColor = baseColor + (1 - baseColor) * blendFraction;
end


function local_format_direction_axis(directionIds, directionLabels)
set(gca, ...
    'XTick', directionIds, ...
    'XTickLabel', directionLabels, ...
    'FontSize', 9);
xlabel('Direction');
xlim([min(directionIds) - 0.5, max(directionIds) + 0.5]);
end
