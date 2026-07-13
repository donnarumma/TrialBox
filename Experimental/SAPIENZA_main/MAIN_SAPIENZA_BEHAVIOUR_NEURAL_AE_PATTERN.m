%% MAIN_SAPIENZA_BEHAVIOUR_NEURAL_AE_PATTERN
% Targeted behaviour-neural checks for the direction-specific |AE| pattern.
%
% This script uses the already-built BehaviourNeural_jointRows.csv table and
% asks whether the direction profile of the behavioural |AE| joint benefit has
% a matching neural profile within each session.

clear;
close all;

analysisPms = struct();
analysisPms.inputRoot = '/TESTS/SAPIENZA/BEHAVIOUR/NeuralCorrelation';
analysisPms.outputRoot = fullfile(analysisPms.inputRoot, 'AEAbsPattern');
analysisPms.visible = 'off';
analysisPms.directionIds = 1:8;
analysisPms.directionLabels = {'D1 forward', 'D2 up-right', 'D3 right', ...
    'D4 down-right', 'D5 backward', 'D6 down-left', 'D7 left', 'D8 up-left'};
analysisPms.minimumSamplesForCorrelation = 5;
analysisPms.signFlipIterations = 20000;
analysisPms.monkeySColor = [0.20 0.40 0.85];
analysisPms.monkeyKColor = [0.20 0.60 0.20];

analysisPms.focusPairs = cell2table({ ...
    'Frontal',  'M7', 'SJointBenefit_AE_abs', 'JointGainSToK',              'S |AE| benefit vs S->K gain'; ...
    'Frontal',  'M8', 'SJointBenefit_AE_abs', 'JointGainSToK',              'S |AE| benefit vs S->K gain'; ...
    'Frontal',  'M4', 'SJointBenefit_AE_abs', 'JointGainSToK',              'S |AE| benefit vs S->K gain'; ...
    'Parietal', 'M4', 'KJointBenefit_AE_abs', 'NeuralJointGainBalance',     'K |AE| benefit vs gain balance'; ...
    'Parietal', 'M7', 'KJointBenefit_AE_abs', 'NeuralJointGainBalance',     'K |AE| benefit vs gain balance'}, ...
    'VariableNames', {'Chamber', 'Model', 'BehaviourReadout', 'NeuralReadout', 'Description'});

if ~exist(analysisPms.outputRoot, 'dir')
    mkdir(analysisPms.outputRoot);
end

plotDir = fullfile(analysisPms.outputRoot, 'plots');
pdfDir = fullfile(plotDir, 'PDF');
figDir = fullfile(plotDir, 'FIG');
if ~exist(pdfDir, 'dir')
    mkdir(pdfDir);
end
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

jointRowsPath = fullfile(analysisPms.inputRoot, 'BehaviourNeural_jointRows.csv');
if ~exist(jointRowsPath, 'file')
    error('%s:MissingInput', mfilename, 'Missing %s', jointRowsPath);
end

jointTable = readtable(jointRowsPath);

[globalTable, directionalTable, sessionPatternTable, sessionSummaryTable] = ...
    local_build_targeted_ae_tables(jointTable, analysisPms);

writetable(globalTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_AEAbs_targetedGlobalCorrelations.csv'));
writetable(directionalTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_AEAbs_targetedDirectionalCorrelations.csv'));
writetable(sessionPatternTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_AEAbs_sessionPatternCorrelations.csv'));
writetable(sessionSummaryTable, fullfile(analysisPms.outputRoot, ...
    'BehaviourNeural_AEAbs_sessionPatternSummary.csv'));

local_plot_targeted_scatter(jointTable, globalTable, analysisPms, pdfDir, figDir);
local_plot_targeted_directional_bars(directionalTable, globalTable, analysisPms, pdfDir, figDir);
local_plot_session_pattern_test(sessionPatternTable, sessionSummaryTable, analysisPms, pdfDir, figDir);

save(fullfile(analysisPms.outputRoot, 'BehaviourNeural_AEAbs_patternAnalysis.mat'), ...
    'analysisPms', 'globalTable', 'directionalTable', ...
    'sessionPatternTable', 'sessionSummaryTable', '-v7.3');

fprintf('Targeted AE_abs behaviour-neural analysis saved in %s\n', analysisPms.outputRoot);


function [globalTable, directionalTable, sessionPatternTable, sessionSummaryTable] = ...
    local_build_targeted_ae_tables(jointTable, analysisPms)

globalRows = {};
directionalRows = {};
sessionRows = {};

for focusPairCounter = 1:height(analysisPms.focusPairs)
    focusPair = analysisPms.focusPairs(focusPairCounter, :);
    pairMask = strcmp(jointTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(jointTable.Model, focusPair.Model{1});
    pairTable = jointTable(pairMask, :);

    behaviourValues = pairTable.(focusPair.BehaviourReadout{1});
    neuralValues = pairTable.(focusPair.NeuralReadout{1});
    [rhoValue, pValue, sampleCount] = local_spearman(behaviourValues, neuralValues, ...
        analysisPms.minimumSamplesForCorrelation);
    globalRows(end + 1, :) = { ...
        focusPair.Chamber{1}, focusPair.Model{1}, focusPair.BehaviourReadout{1}, ...
        focusPair.NeuralReadout{1}, focusPair.Description{1}, sampleCount, ...
        rhoValue, pValue}; %#ok<AGROW>

    for directionCounter = 1:numel(analysisPms.directionIds)
        directionId = analysisPms.directionIds(directionCounter);
        directionTable = pairTable(pairTable.Direction == directionId, :);
        directionBehaviour = directionTable.(focusPair.BehaviourReadout{1});
        directionNeural = directionTable.(focusPair.NeuralReadout{1});
        [directionRho, directionPValue, directionSampleCount] = local_spearman( ...
            directionBehaviour, directionNeural, analysisPms.minimumSamplesForCorrelation);
        directionalRows(end + 1, :) = { ...
            focusPair.Chamber{1}, focusPair.Model{1}, directionId, ...
            analysisPms.directionLabels{directionCounter}, focusPair.BehaviourReadout{1}, ...
            focusPair.NeuralReadout{1}, focusPair.Description{1}, ...
            directionSampleCount, directionRho, directionPValue}; %#ok<AGROW>
    end

    sessionNames = unique(pairTable.Session, 'stable');
    for sessionCounter = 1:numel(sessionNames)
        sessionName = sessionNames{sessionCounter};
        sessionTable = pairTable(strcmp(pairTable.Session, sessionName), :);
        sessionTable = sortrows(sessionTable, 'Direction');
        sessionBehaviour = sessionTable.(focusPair.BehaviourReadout{1});
        sessionNeural = sessionTable.(focusPair.NeuralReadout{1});
        [sessionRho, sessionPValue, sessionSampleCount] = local_spearman( ...
            sessionBehaviour, sessionNeural, analysisPms.minimumSamplesForCorrelation);
        sessionRows(end + 1, :) = { ...
            focusPair.Chamber{1}, sessionName, focusPair.Model{1}, ...
            focusPair.BehaviourReadout{1}, focusPair.NeuralReadout{1}, ...
            focusPair.Description{1}, sessionSampleCount, sessionRho, sessionPValue, ...
            local_direction_profile_string(sessionTable.Direction, sessionBehaviour), ...
            local_direction_profile_string(sessionTable.Direction, sessionNeural)}; %#ok<AGROW>
    end
end

globalTable = cell2table(globalRows, ...
    'VariableNames', {'Chamber', 'Model', 'BehaviourReadout', 'NeuralReadout', ...
    'Description', 'N', 'Rho', 'PValue'});
directionalTable = cell2table(directionalRows, ...
    'VariableNames', {'Chamber', 'Model', 'Direction', 'DirectionLabel', ...
    'BehaviourReadout', 'NeuralReadout', 'Description', 'N', 'Rho', 'PValue'});
sessionPatternTable = cell2table(sessionRows, ...
    'VariableNames', {'Chamber', 'Session', 'Model', 'BehaviourReadout', ...
    'NeuralReadout', 'Description', 'N', 'PatternRho', 'PatternPValue', ...
    'BehaviourProfile', 'NeuralProfile'});
sessionSummaryTable = local_build_session_summary(sessionPatternTable, analysisPms);
end


function sessionSummaryTable = local_build_session_summary(sessionPatternTable, analysisPms)
summaryRows = {};
if isempty(sessionPatternTable)
    sessionSummaryTable = table();
    return;
end

groupKeys = unique(sessionPatternTable(:, ...
    {'Chamber', 'Model', 'BehaviourReadout', 'NeuralReadout', 'Description'}), 'rows');

for groupCounter = 1:height(groupKeys)
    groupMask = strcmp(sessionPatternTable.Chamber, groupKeys.Chamber{groupCounter}) & ...
        strcmp(sessionPatternTable.Model, groupKeys.Model{groupCounter}) & ...
        strcmp(sessionPatternTable.BehaviourReadout, groupKeys.BehaviourReadout{groupCounter}) & ...
        strcmp(sessionPatternTable.NeuralReadout, groupKeys.NeuralReadout{groupCounter});
    groupRhos = sessionPatternTable.PatternRho(groupMask);
    groupRhos = groupRhos(isfinite(groupRhos));
    sessionCount = numel(groupRhos);
    if sessionCount == 0
        meanRho = NaN;
        medianRho = NaN;
        semRho = NaN;
        positiveCount = 0;
        negativeCount = 0;
        signFlipPValue = NaN;
    else
        meanRho = mean(groupRhos, 'omitnan');
        medianRho = median(groupRhos, 'omitnan');
        semRho = std(groupRhos, 'omitnan') ./ sqrt(sessionCount);
        positiveCount = sum(groupRhos > 0);
        negativeCount = sum(groupRhos < 0);
        signFlipPValue = local_signflip_pvalue(groupRhos, analysisPms.signFlipIterations);
    end

    summaryRows(end + 1, :) = { ...
        groupKeys.Chamber{groupCounter}, groupKeys.Model{groupCounter}, ...
        groupKeys.BehaviourReadout{groupCounter}, groupKeys.NeuralReadout{groupCounter}, ...
        groupKeys.Description{groupCounter}, sessionCount, meanRho, medianRho, semRho, ...
        positiveCount, negativeCount, signFlipPValue}; %#ok<AGROW>
end

sessionSummaryTable = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Model', 'BehaviourReadout', 'NeuralReadout', ...
    'Description', 'SessionCount', 'MeanPatternRho', 'MedianPatternRho', ...
    'SemPatternRho', 'PositiveSessionCount', 'NegativeSessionCount', ...
    'SignFlipPValue'});
end


function local_plot_targeted_scatter(jointTable, globalTable, analysisPms, pdfDir, figDir)
figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'Targeted AE_abs behaviour-neural scatter', ...
    'Position', [80 80 1500 900]);
tileLayoutHandle = tiledlayout(2, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');
directionColors = turbo(numel(analysisPms.directionIds));

for focusPairCounter = 1:height(analysisPms.focusPairs)
    focusPair = analysisPms.focusPairs(focusPairCounter, :);
    nexttile;
    pairMask = strcmp(jointTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(jointTable.Model, focusPair.Model{1});
    pairTable = jointTable(pairMask, :);
    behaviourValues = pairTable.(focusPair.BehaviourReadout{1});
    neuralValues = pairTable.(focusPair.NeuralReadout{1});
    directionValues = pairTable.Direction;
    validMask = isfinite(behaviourValues) & isfinite(neuralValues) & isfinite(directionValues);
    scatterColors = directionColors(directionValues(validMask), :);
    scatter(behaviourValues(validMask), neuralValues(validMask), 32, scatterColors, ...
        'filled', 'MarkerFaceAlpha', 0.72, 'MarkerEdgeColor', 'none');
    hold on;
    local_plot_linear_fit(behaviourValues(validMask), neuralValues(validMask));
    grid on;
    colorbarHandle = colorbar;
    colormap(gca, directionColors);
    clim([0.5 8.5]);
    colorbarHandle.Ticks = analysisPms.directionIds;
    colorbarHandle.TickLabels = compose('D%d', analysisPms.directionIds);
    globalMask = strcmp(globalTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(globalTable.Model, focusPair.Model{1}) & ...
        strcmp(globalTable.BehaviourReadout, focusPair.BehaviourReadout{1}) & ...
        strcmp(globalTable.NeuralReadout, focusPair.NeuralReadout{1});
    globalRow = globalTable(globalMask, :);
    title({sprintf('%s %s: %s', ...
        focusPair.Chamber{1}, focusPair.Model{1}, focusPair.Description{1}), ...
        sprintf('rho=%+.2f, p=%s', ...
        globalRow.Rho(1), local_pvalue_text(globalRow.PValue(1)))}, ...
        'Interpreter', 'none');
    xlabel(local_pretty_readout(focusPair.BehaviourReadout{1}), 'Interpreter', 'none');
    ylabel(local_pretty_readout(focusPair.NeuralReadout{1}), 'Interpreter', 'none');
end

title(tileLayoutHandle, ...
    'Targeted |AE| behaviour-neural associations across sessions and directions', ...
    'Interpreter', 'none');
figureStem = 'BehaviourNeural_AEAbs_TargetedScatters';
exportgraphics(figureHandle, fullfile(pdfDir, [figureStem '.pdf']), ...
    'ContentType', 'vector');
savefig(figureHandle, fullfile(figDir, [figureStem '.fig']));
close(figureHandle);
end


function local_plot_targeted_directional_bars(directionalTable, globalTable, analysisPms, pdfDir, figDir)
figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'Targeted AE_abs directional correlations', ...
    'Position', [80 80 1500 900]);
tileLayoutHandle = tiledlayout(2, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for focusPairCounter = 1:height(analysisPms.focusPairs)
    focusPair = analysisPms.focusPairs(focusPairCounter, :);
    nexttile;
    directionMask = strcmp(directionalTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(directionalTable.Model, focusPair.Model{1}) & ...
        strcmp(directionalTable.BehaviourReadout, focusPair.BehaviourReadout{1}) & ...
        strcmp(directionalTable.NeuralReadout, focusPair.NeuralReadout{1});
    pairDirectionalTable = directionalTable(directionMask, :);
    pairDirectionalTable = sortrows(pairDirectionalTable, 'Direction');

    barColor = local_pair_color(focusPair.BehaviourReadout{1}, analysisPms);
    barHandle = bar(pairDirectionalTable.Direction, pairDirectionalTable.Rho, ...
        'FaceColor', 'flat', 'EdgeColor', 'none');
    for directionCounter = 1:height(pairDirectionalTable)
        if pairDirectionalTable.Rho(directionCounter) >= 0
            barHandle.CData(directionCounter, :) = barColor;
        else
            barHandle.CData(directionCounter, :) = local_lighten_color(barColor, 0.55);
        end
    end
    hold on;
    yline(0, 'k-');
    globalMask = strcmp(globalTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(globalTable.Model, focusPair.Model{1}) & ...
        strcmp(globalTable.BehaviourReadout, focusPair.BehaviourReadout{1}) & ...
        strcmp(globalTable.NeuralReadout, focusPair.NeuralReadout{1});
    globalRow = globalTable(globalMask, :);
    yline(globalRow.Rho(1), '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
    grid on;
    xticks(analysisPms.directionIds);
    xticklabels(compose('D%d', analysisPms.directionIds));
    ylim(local_symmetric_limits([pairDirectionalTable.Rho; globalRow.Rho(1)]));
    title(sprintf('%s %s: %s', focusPair.Chamber{1}, focusPair.Model{1}, ...
        focusPair.Description{1}), 'Interpreter', 'none');
    ylabel('Direction-wise Spearman \rho');
    xlabel('Direction');
    for directionCounter = 1:height(pairDirectionalTable)
        rhoValue = pairDirectionalTable.Rho(directionCounter);
        if isfinite(rhoValue)
            local_add_bar_label(pairDirectionalTable.Direction(directionCounter), ...
                rhoValue, local_pvalue_text(pairDirectionalTable.PValue(directionCounter)));
        end
    end
end

title(tileLayoutHandle, ...
    'Targeted direction-wise correlations for |AE| behaviour-neural links', ...
    'Interpreter', 'none');
figureStem = 'BehaviourNeural_AEAbs_TargetedDirectionalBars';
exportgraphics(figureHandle, fullfile(pdfDir, [figureStem '.pdf']), ...
    'ContentType', 'vector');
savefig(figureHandle, fullfile(figDir, [figureStem '.fig']));
close(figureHandle);
end


function local_plot_session_pattern_test(sessionPatternTable, sessionSummaryTable, analysisPms, pdfDir, figDir)
figureHandle = figure( ...
    'Color', 'w', ...
    'Visible', analysisPms.visible, ...
    'Name', 'Targeted AE_abs session pattern test', ...
    'Position', [80 80 1450 760]);
hold on;

pairCount = height(analysisPms.focusPairs);
for focusPairCounter = 1:pairCount
    focusPair = analysisPms.focusPairs(focusPairCounter, :);
    pairMask = strcmp(sessionPatternTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(sessionPatternTable.Model, focusPair.Model{1}) & ...
        strcmp(sessionPatternTable.BehaviourReadout, focusPair.BehaviourReadout{1}) & ...
        strcmp(sessionPatternTable.NeuralReadout, focusPair.NeuralReadout{1});
    pairTable = sessionPatternTable(pairMask, :);
    pairTable = sortrows(pairTable, 'Session');
    pairColor = local_pair_color(focusPair.BehaviourReadout{1}, analysisPms);
    jitterOffsets = linspace(-0.22, 0.22, max(height(pairTable), 2))';
    if height(pairTable) == 1
        jitterOffsets = 0;
    else
        jitterOffsets = jitterOffsets(randperm(numel(jitterOffsets)));
    end
    scatter(focusPairCounter + jitterOffsets(1:height(pairTable)), ...
        pairTable.PatternRho, 44, pairColor, 'filled', ...
        'MarkerFaceAlpha', 0.72, 'MarkerEdgeColor', 'w');

    summaryMask = strcmp(sessionSummaryTable.Chamber, focusPair.Chamber{1}) & ...
        strcmp(sessionSummaryTable.Model, focusPair.Model{1}) & ...
        strcmp(sessionSummaryTable.BehaviourReadout, focusPair.BehaviourReadout{1}) & ...
        strcmp(sessionSummaryTable.NeuralReadout, focusPair.NeuralReadout{1});
    summaryRow = sessionSummaryTable(summaryMask, :);
    errorbar(focusPairCounter, summaryRow.MeanPatternRho(1), ...
        summaryRow.SemPatternRho(1), 'ks', ...
        'MarkerFaceColor', 'w', 'LineWidth', 1.4, 'CapSize', 10);
    text(focusPairCounter, -0.92, sprintf('p=%s', ...
        local_pvalue_text(summaryRow.SignFlipPValue(1))), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

yline(0, 'k-');
grid on;
xlim([0.5 pairCount + 0.5]);
ylim([-1 1]);
xticks(1:pairCount);
xticklabels(local_focus_pair_labels(analysisPms.focusPairs));
xtickangle(25);
ylabel('Within-session Spearman \rho across directions');
title({'Session-wise pattern test', ...
    'Each point correlates the 8-direction behavioural |AE| profile with the matching neural profile within one session'}, ...
    'Interpreter', 'none');
figureStem = 'BehaviourNeural_AEAbs_SessionPatternTest';
exportgraphics(figureHandle, fullfile(pdfDir, [figureStem '.pdf']), ...
    'ContentType', 'vector');
savefig(figureHandle, fullfile(figDir, [figureStem '.fig']));
close(figureHandle);
end


function [rhoValue, pValue, sampleCount] = local_spearman(firstValues, secondValues, minimumSamples)
validMask = isfinite(firstValues) & isfinite(secondValues);
sampleCount = sum(validMask);
if sampleCount < minimumSamples
    rhoValue = NaN;
    pValue = NaN;
    return;
end
[rhoValue, pValue] = corr(firstValues(validMask), secondValues(validMask), ...
    'Type', 'Spearman', 'Rows', 'complete');
end


function pValue = local_signflip_pvalue(values, iterationCount)
values = values(isfinite(values));
sampleCount = numel(values);
if sampleCount < 2
    pValue = NaN;
    return;
end

observedMean = mean(values);
if observedMean == 0
    pValue = 1;
    return;
end

rng(7741, 'twister');
randomSigns = (rand(iterationCount, sampleCount) > 0.5) .* 2 - 1;
randomMeans = mean(randomSigns .* values(:)', 2);
pValue = mean(abs(randomMeans) >= abs(observedMean));
pValue = max(pValue, 1 ./ iterationCount);
end


function profileText = local_direction_profile_string(directionValues, readoutValues)
profileParts = strings(numel(directionValues), 1);
for directionCounter = 1:numel(directionValues)
    profileParts(directionCounter) = sprintf('D%d=%+.3f', ...
        directionValues(directionCounter), readoutValues(directionCounter));
end
profileText = strjoin(profileParts, ' ');
end


function local_plot_linear_fit(behaviourValues, neuralValues)
validMask = isfinite(behaviourValues) & isfinite(neuralValues);
if sum(validMask) < 3
    return;
end
fitCoefficients = polyfit(behaviourValues(validMask), neuralValues(validMask), 1);
xLimits = xlim;
fitX = linspace(xLimits(1), xLimits(2), 100);
fitY = polyval(fitCoefficients, fitX);
plot(fitX, fitY, 'k-', 'LineWidth', 1.1);
end


function readoutLabel = local_pretty_readout(readoutName)
readoutLabel = strrep(readoutName, '_', ' ');
readoutLabel = strrep(readoutLabel, 'AE abs', '|AE|');
readoutLabel = strrep(readoutLabel, 'SToK', 'S->K');
readoutLabel = strrep(readoutLabel, 'KToS', 'K->S');
end


function pValueText = local_pvalue_text(pValue)
if ~isfinite(pValue)
    pValueText = 'n/a';
elseif pValue < 0.001
    pValueText = sprintf('%.1e', pValue);
else
    pValueText = sprintf('%.3f', pValue);
end
end


function axisLimits = local_symmetric_limits(values)
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    axisLimits = [-1 1];
    return;
end
maxAbsValue = max(abs(finiteValues));
maxAbsValue = max(maxAbsValue, 0.1);
axisLimits = [-1.12 .* maxAbsValue, 1.12 .* maxAbsValue];
end


function local_add_bar_label(directionId, rhoValue, pValueText)
if rhoValue >= 0
    verticalAlignment = 'bottom';
    yOffset = 0.03;
else
    verticalAlignment = 'top';
    yOffset = -0.03;
end
text(directionId, rhoValue + yOffset, pValueText, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', verticalAlignment, ...
    'FontSize', 8, ...
    'Rotation', 90);
end


function labels = local_focus_pair_labels(focusPairs)
labels = strings(height(focusPairs), 1);
for focusPairCounter = 1:height(focusPairs)
    if startsWith(focusPairs.BehaviourReadout{focusPairCounter}, 'S')
        receiverTag = 'S';
    else
        receiverTag = 'K';
    end
    labels(focusPairCounter) = sprintf('%s %s %s', ...
        focusPairs.Chamber{focusPairCounter}, focusPairs.Model{focusPairCounter}, receiverTag);
end
end


function pairColor = local_pair_color(behaviourReadout, analysisPms)
if startsWith(behaviourReadout, 'S')
    pairColor = analysisPms.monkeySColor;
else
    pairColor = analysisPms.monkeyKColor;
end
end


function lightColor = local_lighten_color(baseColor, blendFraction)
lightColor = baseColor + (1 - baseColor) .* blendFraction;
end
