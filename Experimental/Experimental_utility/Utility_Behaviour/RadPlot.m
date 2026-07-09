%function RadPlot(dataS,dataK,par)
function RadPlot(dataS,dataK,par)

legendNames = par.legendNames;
ind_correctDir = [3 2 1 8 7 6 5 4]; % Rimappatura spaziale corretta

if contains(par.InField, {'AE', 'angular','Hasson'}, 'IgnoreCase', true)
    S_Rad = abs(dataS(1,ind_correctDir)); 
    K_Rad = abs(dataK(1,ind_correctDir));
else
    S_Rad = dataS(1,ind_correctDir);
    K_Rad = dataK(1,ind_correctDir);
end

if strcmp(par.mod,'Zero')
    S_Rad(isnan(S_Rad)) = 0;
    K_Rad(isnan(K_Rad)) = 0;
elseif strcmp(par.mod,'Interp')
    if any(isnan(S_Rad))
        S_Rad = fillmissing(S_Rad, 'linear', 'EndValues', 'nearest');
    end
    if any(isnan(K_Rad))
        K_Rad = fillmissing(K_Rad, 'linear', 'EndValues', 'nearest');
    end
end

max_v = par.maxValue;
if isfield(par, 'minValue')
    min_v = par.minValue;
else
    min_v = 0;
end

if min_v < 0
    min_v = 0;
end

figure('Name', [par.InField,'S vs K'], 'Color', 'w');
theta = linspace(0, 2*pi, 9);
nexttile;

% === FONT FLESSIBILI (STESSO DI RadMultiplePlot) ===
fs_ticks_default = 10;
fs_legend_default = 10;
if isfield(par, 'font_ticks'),  fs_ticks = par.font_ticks; else, fs_ticks = fs_ticks_default; end
if isfield(par, 'font_legend'), fs_legend = par.font_legend; else, fs_legend = fs_legend_default; end
if isfield(par, 'font_title'),  fs_title_mult = par.font_title / fs_ticks; else, fs_title_mult = 1.2; end
% === FINE FONT ===

s_p = [S_Rad(1,:), S_Rad(1,1)];
k_p = [K_Rad(1,:), K_Rad(1,1)];

colors = {par.color1, par.color2};
styles = {par.lineStyle1, par.lineStyle2};

hS = polarplot(theta, s_p, styles{1}, 'LineWidth', 2, ...
    'Color', colors{1}, 'HandleVisibility', 'on'); hold on;
polarplot(theta, s_p, 'o', 'MarkerFaceColor', colors{1}, ...
    'MarkerEdgeColor', colors{1}, 'HandleVisibility', 'off', 'MarkerSize', 4);

hK = polarplot(theta, k_p, styles{2}, 'LineWidth', 2, ...
    'Color', colors{2}, 'HandleVisibility', 'on'); hold on;
polarplot(theta, k_p, 'o', 'MarkerFaceColor', colors{2}, ...
    'MarkerEdgeColor', colors{2}, 'HandleVisibility', 'off', 'MarkerSize', 4);

% Stelle p-value
if isfield(par, 'p_vals')
    p_vals_curr = par.p_vals(1, ind_correctDir);
    r_star = max_v * 0.95;
    for d = 1:8
        if p_vals_curr(d) <= 0.001           
            text(theta(d), r_star, '***', 'Color', 'r', 'FontSize', 18, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'HandleVisibility', 'off');
        elseif p_vals_curr(d) <= 0.01            
            text(theta(d), r_star, '**', 'Color', 'r', 'FontSize', 18, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'HandleVisibility', 'off');
        elseif p_vals_curr(d) <= 0.05
            text(theta(d), r_star, '*', 'Color', 'r', 'FontSize', 18, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'HandleVisibility', 'off');
        end
    end
end

% Titoli (corretto sprintf)
colorStr1 = sprintf('%.2g,%.2g,%.2g', par.color1);
titlerow1 = ['\color[rgb]{', colorStr1, '}', par.cond_names1{1}];
colorStr2 = sprintf('%.2g,%.2g,%.2g', par.color2);
titlerow2 = ['\color[rgb]{', colorStr2, '}', par.cond_names2{1}];
pax = gca;
title(pax, {titlerow1, titlerow2}, 'Interpreter', 'tex');
if isfield(par, 'unitMeasure')
    pax.RAxis.Label.String = par.unitMeasure;
    pax.RAxis.Label.Interpreter = 'tex';
end
% Assi + FONT applicati
pax.TitleHorizontalAlignment = 'center';
pax.FontSize = fs_ticks;        % ← Override applicato
pax.FontWeight = 'bold';
pax.TitleFontSizeMultiplier = fs_title_mult;  % ← Override applicato
pax.ThetaTick = linspace(0, 315, 8);
pax.ThetaTickLabel = {'D3','D2','D1','D8','D7','D6','D5','D4'};
pax.RLim = [min_v max_v];       % ← Zoom!

% Leggenda con override
legendNames_plot = cellfun(@(c) strjoin(c, '\n '), legendNames, 'UniformOutput', false);
lgd = legend([hS, hK], legendNames_plot, 'Orientation', 'horizontal', 'Interpreter', 'tex'); 
lgd.FontSize = fs_legend;       % ← Override applicato
lgd.Layout.Tile = 'south';

set(gcf, 'WindowState', 'maximized');
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])

% Salvataggio
pMaxDistDiff_path = [par.save_dir,filesep];
par.savePlotEpsPdfMat.dir_png = fullfile(pMaxDistDiff_path, 'PNGs');
par.savePlotEpsPdfMat.dir_pdf = fullfile(pMaxDistDiff_path, 'PDFs');
par.savePlotEpsPdfMat.dir_mat = fullfile(pMaxDistDiff_path, 'MATfiles');
name_fig =['PolarPlot_',par.InField,'_',par.figureName];
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
close(gcf)
end
