%% MAIN_SAPIENZA_BEHAVIOUR_RT_AE_PATTERN
% Direction-resolved RT and |AE| tradeoff in source-active versus joint action.
%
% Positive benefit means joint action improves the metric relative to the
% source-active baseline. For RT and |AE| this means lower values in joint.

clear;
close all;

analysisPms = struct();
analysisPms.inputRoot = '/TESTS/SAPIENZA/BEHAVIOUR/DirectionalPlots';
analysisPms.outputRoot = fullfile(analysisPms.inputRoot, 'RTTradeoff');
analysisPms.summaryPath = fullfile(analysisPms.inputRoot, 'BehaviourDirectional_summary.csv');
analysisPms.visible = 'off';
analysisPms.chamberNames = {'Frontal', 'Parietal'};
analysisPms.directionIds = 1:8;
analysisPms.directionLabels = {'D1 forward', 'D2 up-right', 'D3 right', ...
    'D4 down-right', 'D5 backward', 'D6 down-left', 'D7 left', 'D8 up-left'};
analysisPms.directionAnglesDegrees = 0:45:315;
analysisPms.rtPolarMinRadius = 150;
analysisPms.monkeySColor = [0.20 0.40 0.85];
analysisPms.monkeyKColor = [0.20 0.60 0.20];

if ~exist(analysisPms.summaryPath, 'file')
    error('%s:MissingInput', mfilename, 'Missing %s', analysisPms.summaryPath);
end

if ~exist(analysisPms.outputRoot, 'dir')
    mkdir(analysisPms.outputRoot);
end

pdfDir = fullfile(analysisPms.outputRoot, 'plots', 'PDF');
figDir = fullfile(analysisPms.outputRoot, 'plots', 'FIG');
if ~exist(pdfDir, 'dir')
    mkdir(pdfDir);
end
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

summaryTable = readtable(analysisPms.summaryPath);
tradeoffTable = local_build_rt_ae_tradeoff_table(summaryTable, analysisPms);
writetable(tradeoffTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourDirectional_RTAE_Tradeoff.csv'));

local_plot_rt_source_active_vs_joint(summaryTable, analysisPms, pdfDir, figDir);
local_plot_rt_polar_pattern(tradeoffTable, analysisPms, pdfDir, figDir);
local_plot_rt_ae_tradeoff_quadrants(tradeoffTable, analysisPms, pdfDir, figDir);
local_plot_rt_ae_tradeoff_bars(tradeoffTable, analysisPms, pdfDir, figDir);

save(fullfile(analysisPms.outputRoot, 'BehaviourDirectional_RTAE_Tradeoff.mat'), ...
    'analysisPms', 'tradeoffTable', '-v7.3');

fprintf('RT/AE directional tradeoff analysis saved in %s\n', analysisPms.outputRoot);


function tradeoffTable = local_build_rt_ae_tradeoff_table(summaryTable, analysisPms)
tradeoffRows = {};
rtTable = summaryTable(strcmp(summaryTable.Metric, 'RT'), :);
aeTable = summaryTable(strcmp(summaryTable.Metric, 'AE_abs'), :);

for chamberCounter = 1:numel(analysisPms.chamberNames)
    chamberName = analysisPms.chamberNames{chamberCounter};
    for directionCounter = 1:numel(analysisPms.directionIds)
        directionId = analysisPms.directionIds(directionCounter);
        rtRow = rtTable(strcmp(rtTable.Chamber, chamberName) & rtTable.Direction == directionId, :);
        aeRow = aeTable(strcmp(aeTable.Chamber, chamberName) & aeTable.Direction == directionId, :);
        if isempty(rtRow) || isempty(aeRow)
            continue;
        end

        tradeoffRows(end + 1, :) = local_tradeoff_row( ...
            chamberName, 'S', directionId, analysisPms.directionLabels{directionCounter}, ...
            rtRow.MeanSAct(1), rtRow.SemSAct(1), rtRow.MeanSJoint(1), rtRow.SemSJoint(1), ...
            rtRow.MeanSJointBenefit(1), rtRow.SemSJointBenefit(1), ...
            aeRow.MeanSAct(1), aeRow.SemSAct(1), aeRow.MeanSJoint(1), aeRow.SemSJoint(1), ...
            aeRow.MeanSJointBenefit(1), aeRow.SemSJointBenefit(1)); %#ok<AGROW>

        tradeoffRows(end + 1, :) = local_tradeoff_row( ...
            chamberName, 'K', directionId, analysisPms.directionLabels{directionCounter}, ...
            rtRow.MeanKAct(1), rtRow.SemKAct(1), rtRow.MeanKJoint(1), rtRow.SemKJoint(1), ...
            rtRow.MeanKJointBenefit(1), rtRow.SemKJointBenefit(1), ...
            aeRow.MeanKAct(1), aeRow.SemKAct(1), aeRow.MeanKJoint(1), aeRow.SemKJoint(1), ...
            aeRow.MeanKJointBenefit(1), aeRow.SemKJointBenefit(1)); %#ok<AGROW>
    end
end

tradeoffTable = cell2table(tradeoffRows, ...
    'VariableNames', {'Chamber', 'Monkey', 'Direction', 'DirectionLabel', ...
    'RTSourceActive', 'SemRTSourceActive', 'RTJoint', 'SemRTJoint', ...
    'RTBenefit', 'SemRTBenefit', 'AEAbsSourceActive', 'SemAEAbsSourceActive', ...
    'AEAbsJoint', 'SemAEAbsJoint', 'AEAbsBenefit', 'SemAEAbsBenefit', ...
    'TradeoffClass'});
end


function rowValues = local_tradeoff_row(chamberName, monkeyName, directionId, directionLabel, ...
    rtSourceActive, semRtSourceActive, rtJoint, semRtJoint, rtBenefit, semRtBenefit, ...
    aeSourceActive, semAeSourceActive, aeJoint, semAeJoint, aeBenefit, semAeBenefit)

tradeoffClass = local_tradeoff_class(rtBenefit, aeBenefit);
rowValues = {chamberName, monkeyName, directionId, directionLabel, ...
    rtSourceActive, semRtSourceActive, rtJoint, semRtJoint, rtBenefit, semRtBenefit, ...
    aeSourceActive, semAeSourceActive, aeJoint, semAeJoint, aeBenefit, semAeBenefit, ...
    tradeoffClass};
end


function tradeoffClass = local_tradeoff_class(rtBenefit, aeBenefit)
if rtBenefit >= 0 && aeBenefit >= 0
    tradeoffClass = 'faster and more accurate';
elseif rtBenefit < 0 && aeBenefit >= 0
    tradeoffClass = 'slower but more accurate';
elseif rtBenefit >= 0 && aeBenefit < 0
    tradeoffClass = 'faster but less accurate';
else
    tradeoffClass = 'slower and less accurate';
end
end


function local_plot_rt_source_active_vs_joint(summaryTable, analysisPms, pdfDir, figDir)
rtTable = summaryTable(strcmp(summaryTable.Metric, 'RT'), :);
directionIds = analysisPms.directionIds;
directionLabels = compose('D%d', directionIds);
lightSColor = local_lighten_color(analysisPms.monkeySColor, 0.58);
lightKColor = local_lighten_color(analysisPms.monkeyKColor, 0.58);

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'RT source-active versus joint by direction', ...
    'Position', [80 80 1460 760]);
tiledlayout(numel(analysisPms.chamberNames), 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for chamberCounter = 1:numel(analysisPms.chamberNames)
    chamberName = analysisPms.chamberNames{chamberCounter};
    chamberTable = rtTable(strcmp(rtTable.Chamber, chamberName), :);
    chamberTable = sortrows(chamberTable, 'Direction');

    nexttile;
    local_plot_condition_pair(directionIds, ...
        chamberTable.MeanSAct, chamberTable.SemSAct, ...
        chamberTable.MeanSJoint, chamberTable.SemSJoint, ...
        lightSColor, analysisPms.monkeySColor, {'S ActS,ObsK', 'S Joint'});
    title(sprintf('%s S: source-active vs joint', chamberName), 'Interpreter', 'none');
    ylabel('RT (lower is faster)');
    local_format_direction_axis(directionIds, directionLabels);

    nexttile;
    local_plot_condition_pair(directionIds, ...
        chamberTable.MeanKAct, chamberTable.SemKAct, ...
        chamberTable.MeanKJoint, chamberTable.SemKJoint, ...
        lightKColor, analysisPms.monkeyKColor, {'K ObsS,ActK', 'K Joint'});
    title(sprintf('%s K: source-active vs joint', chamberName), 'Interpreter', 'none');
    ylabel('RT (lower is faster)');
    local_format_direction_axis(directionIds, directionLabels);

    nexttile;
    local_plot_two_series(directionIds, ...
        chamberTable.MeanSJointBenefit, chamberTable.SemSJointBenefit, ...
        chamberTable.MeanKJointBenefit, chamberTable.SemKJointBenefit, ...
        analysisPms.monkeySColor, analysisPms.monkeyKColor);
    yline(0, 'k-');
    title(sprintf('%s: joint benefit in RT', chamberName), 'Interpreter', 'none');
    ylabel('Benefit = source-active - joint');
    if chamberCounter == 1
        legend({'S benefit', 'K benefit'}, 'Location', 'best');
    end
    local_format_direction_axis(directionIds, directionLabels);
end

sgtitle('Reaction time by condition and direction', 'Interpreter', 'none');
figureStem = 'BehaviourDirectional_RT_SourceActiveVsJoint';
exportgraphics(figureHandle, fullfile(pdfDir, [figureStem '.pdf']), 'ContentType', 'vector');
savefig(figureHandle, fullfile(figDir, [figureStem '.fig']));
close(figureHandle);
end


function local_plot_rt_ae_tradeoff_quadrants(tradeoffTable, analysisPms, pdfDir, figDir)
figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'RT and AE tradeoff quadrants', ...
    'Position', [80 80 1240 980]);
tiledlayout(numel(analysisPms.chamberNames), 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

xLimits = local_symmetric_limits(tradeoffTable.RTBenefit);
yLimits = local_symmetric_limits(tradeoffTable.AEAbsBenefit);

for chamberCounter = 1:numel(analysisPms.chamberNames)
    chamberName = analysisPms.chamberNames{chamberCounter};
    for monkeyCounter = 1:2
        if monkeyCounter == 1
            monkeyName = 'S';
            monkeyColor = analysisPms.monkeySColor;
        else
            monkeyName = 'K';
            monkeyColor = analysisPms.monkeyKColor;
        end

        nexttile;
        rowMask = strcmp(tradeoffTable.Chamber, chamberName) & ...
            strcmp(tradeoffTable.Monkey, monkeyName);
        plotTable = sortrows(tradeoffTable(rowMask, :), 'Direction');
        scatter(plotTable.RTBenefit, plotTable.AEAbsBenefit, 110, monkeyColor, ...
            'filled', 'MarkerFaceAlpha', 0.78, 'MarkerEdgeColor', 'w');
        hold on;
        plot(plotTable.RTBenefit, plotTable.AEAbsBenefit, '-', ...
            'Color', local_lighten_color(monkeyColor, 0.25), 'LineWidth', 1.2);
        xline(0, 'k-');
        yline(0, 'k-');
        grid on;
        xlim(xLimits);
        ylim(yLimits);
        xlabel('RT benefit: source-active - joint (positive = faster in joint)');
        ylabel('|AE| benefit: source-active - joint (positive = more accurate in joint)');
        title(sprintf('%s %s: speed-accuracy tradeoff', chamberName, monkeyName), ...
            'Interpreter', 'none');

        local_add_quadrant_labels(xLimits, yLimits);
        for rowCounter = 1:height(plotTable)
            text(plotTable.RTBenefit(rowCounter), plotTable.AEAbsBenefit(rowCounter), ...
                sprintf(' D%d', plotTable.Direction(rowCounter)), ...
                'FontWeight', 'bold', 'FontSize', 9, ...
                'Color', [0.05 0.05 0.05]);
        end
    end
end

sgtitle('RT-|AE| joint-action tradeoff by direction', 'Interpreter', 'none');
figureStem = 'BehaviourDirectional_RTAE_TradeoffQuadrants';
exportgraphics(figureHandle, fullfile(pdfDir, [figureStem '.pdf']), 'ContentType', 'vector');
savefig(figureHandle, fullfile(figDir, [figureStem '.fig']));
close(figureHandle);
end


function local_plot_rt_polar_pattern(tradeoffTable, analysisPms, pdfDir, figDir)
polarTable = sortrows(tradeoffTable, {'Chamber', 'Monkey', 'Direction'});
writetable(polarTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourDirectional_RT_PolarPattern.csv'));

thetaRadians = deg2rad([analysisPms.directionAnglesDegrees(:); 360]);
thetaTickLabels = compose('D%d', analysisPms.directionIds);
lightSColor = local_lighten_color(analysisPms.monkeySColor, 0.58);
lightKColor = local_lighten_color(analysisPms.monkeyKColor, 0.58);
maxRadius = max([polarTable.RTSourceActive; polarTable.RTJoint], [], 'omitnan');
maxRadius = max(maxRadius, analysisPms.rtPolarMinRadius + 50);
maxRadius = ceil(maxRadius / 50) * 50;

figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'RT polar behavioural pattern', ...
    'Position', [80 80 1320 1040]);
tileLayoutHandle = tiledlayout(numel(analysisPms.chamberNames), 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for chamberCounter = 1:numel(analysisPms.chamberNames)
    chamberName = analysisPms.chamberNames{chamberCounter};

    sPolarAxes = polaraxes(tileLayoutHandle);
    sPolarAxes.Layout.Tile = (chamberCounter - 1) * 2 + 1;
    sMask = strcmp(polarTable.Chamber, chamberName) & strcmp(polarTable.Monkey, 'S');
    sTable = sortrows(polarTable(sMask, :), 'Direction');
    local_plot_rt_polar_axes(sPolarAxes, thetaRadians, ...
        sTable.RTSourceActive, sTable.RTJoint, ...
        lightSColor, analysisPms.monkeySColor, ...
        analysisPms.rtPolarMinRadius, maxRadius, thetaTickLabels, ...
        sprintf('%s S: RT compass', chamberName), {'S ActS,ObsK', 'S Joint'});

    kPolarAxes = polaraxes(tileLayoutHandle);
    kPolarAxes.Layout.Tile = (chamberCounter - 1) * 2 + 2;
    kMask = strcmp(polarTable.Chamber, chamberName) & strcmp(polarTable.Monkey, 'K');
    kTable = sortrows(polarTable(kMask, :), 'Direction');
    local_plot_rt_polar_axes(kPolarAxes, thetaRadians, ...
        kTable.RTSourceActive, kTable.RTJoint, ...
        lightKColor, analysisPms.monkeyKColor, ...
        analysisPms.rtPolarMinRadius, maxRadius, thetaTickLabels, ...
        sprintf('%s K: RT compass', chamberName), {'K ObsS,ActK', 'K Joint'});
end

title(tileLayoutHandle, ...
    'Polar view of reaction time by movement direction', ...
    'Interpreter', 'none');
figureStem = 'BehaviourDirectional_RT_PolarPattern';
exportgraphics(figureHandle, fullfile(pdfDir, [figureStem '.pdf']), 'ContentType', 'vector');
savefig(figureHandle, fullfile(figDir, [figureStem '.fig']));
close(figureHandle);
end


function local_plot_rt_polar_axes(polarAxesHandle, thetaRadians, ...
    sourceActiveValues, jointValues, sourceActiveColor, jointColor, ...
    minRadius, maxRadius, thetaTickLabels, plotTitle, legendLabels)
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
polarAxesHandle.RLim = [minRadius maxRadius];
polarAxesHandle.RAxis.Label.String = 'RT';
polarAxesHandle.FontSize = 10;
title(polarAxesHandle, plotTitle, ...
    'Interpreter', 'none', ...
    'FontWeight', 'bold');
legend(polarAxesHandle, legendLabels, ...
    'Location', 'southoutside', ...
    'Interpreter', 'none');
end


function local_plot_rt_ae_tradeoff_bars(tradeoffTable, analysisPms, pdfDir, figDir)
figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'RT and AE normalized benefit bars', ...
    'Position', [80 80 1460 820]);
tiledlayout(numel(analysisPms.chamberNames), 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for chamberCounter = 1:numel(analysisPms.chamberNames)
    chamberName = analysisPms.chamberNames{chamberCounter};
    for monkeyCounter = 1:2
        if monkeyCounter == 1
            monkeyName = 'S';
            monkeyColor = analysisPms.monkeySColor;
        else
            monkeyName = 'K';
            monkeyColor = analysisPms.monkeyKColor;
        end

        nexttile;
        rowMask = strcmp(tradeoffTable.Chamber, chamberName) & ...
            strcmp(tradeoffTable.Monkey, monkeyName);
        plotTable = sortrows(tradeoffTable(rowMask, :), 'Direction');
        normalizedRtBenefit = local_zscore(plotTable.RTBenefit);
        normalizedAeBenefit = local_zscore(plotTable.AEAbsBenefit);

        bar(plotTable.Direction - 0.17, normalizedRtBenefit, 0.32, ...
            'FaceColor', local_lighten_color(monkeyColor, 0.45), 'EdgeColor', 'none');
        hold on;
        bar(plotTable.Direction + 0.17, normalizedAeBenefit, 0.32, ...
            'FaceColor', monkeyColor, 'EdgeColor', 'none');
        yline(0, 'k-');
        grid on;
        xticks(analysisPms.directionIds);
        xticklabels(compose('D%d', analysisPms.directionIds));
        ylabel('Direction-normalized benefit (z)');
        xlabel('Direction');
        title(sprintf('%s %s: normalized RT and |AE| benefits', chamberName, monkeyName), ...
            'Interpreter', 'none');
        if chamberCounter == 1 && monkeyCounter == 1
            legend({'RT benefit', '|AE| benefit'}, 'Location', 'best');
        end
    end
end

sgtitle('Normalized joint benefit profiles for RT and |AE|', 'Interpreter', 'none');
figureStem = 'BehaviourDirectional_RTAE_TradeoffBars';
exportgraphics(figureHandle, fullfile(pdfDir, [figureStem '.pdf']), 'ContentType', 'vector');
savefig(figureHandle, fullfile(figDir, [figureStem '.fig']));
close(figureHandle);
end


function local_plot_condition_pair(directionIds, meanFirst, semFirst, meanSecond, semSecond, ...
    firstColor, secondColor, legendLabels)
hold on;
errorbar(directionIds, meanFirst, semFirst, '-o', ...
    'Color', firstColor, 'MarkerFaceColor', firstColor, ...
    'LineWidth', 1.5, 'CapSize', 4);
errorbar(directionIds, meanSecond, semSecond, '-o', ...
    'Color', secondColor, 'MarkerFaceColor', secondColor, ...
    'LineWidth', 1.5, 'CapSize', 4);
grid on;
legend(legendLabels, 'Location', 'best');
end


function local_plot_two_series(directionIds, meanFirst, semFirst, meanSecond, semSecond, ...
    firstColor, secondColor)
hold on;
errorbar(directionIds, meanFirst, semFirst, '-o', ...
    'Color', firstColor, 'MarkerFaceColor', firstColor, ...
    'LineWidth', 1.5, 'CapSize', 4);
errorbar(directionIds, meanSecond, semSecond, '-o', ...
    'Color', secondColor, 'MarkerFaceColor', secondColor, ...
    'LineWidth', 1.5, 'CapSize', 4);
grid on;
end


function local_format_direction_axis(directionIds, directionLabels)
xticks(directionIds);
xticklabels(directionLabels);
xlim([min(directionIds) - 0.35, max(directionIds) + 0.35]);
xlabel('Direction');
end


function axisLimits = local_symmetric_limits(values)
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    axisLimits = [-1 1];
    return;
end
maxAbsValue = max(abs(finiteValues));
maxAbsValue = max(maxAbsValue, 1);
axisLimits = [-1.18 * maxAbsValue, 1.18 * maxAbsValue];
end


function local_add_quadrant_labels(xLimits, yLimits)
text(xLimits(2) * 0.96, yLimits(2) * 0.92, 'faster + accurate', ...
    'HorizontalAlignment', 'right', 'Color', [0.15 0.15 0.15], 'FontSize', 8);
text(xLimits(1) * 0.96, yLimits(2) * 0.92, 'slower + accurate', ...
    'HorizontalAlignment', 'left', 'Color', [0.15 0.15 0.15], 'FontSize', 8);
text(xLimits(2) * 0.96, yLimits(1) * 0.92, 'faster + less accurate', ...
    'HorizontalAlignment', 'right', 'Color', [0.15 0.15 0.15], 'FontSize', 8);
text(xLimits(1) * 0.96, yLimits(1) * 0.92, 'slower + less accurate', ...
    'HorizontalAlignment', 'left', 'Color', [0.15 0.15 0.15], 'FontSize', 8);
end


function zValues = local_zscore(values)
finiteValues = values(isfinite(values));
if numel(finiteValues) < 2 || std(finiteValues) == 0
    zValues = zeros(size(values));
    zValues(~isfinite(values)) = NaN;
    return;
end
zValues = (values - mean(finiteValues)) ./ std(finiteValues);
end


function lightColor = local_lighten_color(baseColor, blendFraction)
lightColor = baseColor + (1 - baseColor) .* blendFraction;
end
