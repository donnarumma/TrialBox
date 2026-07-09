function RadPlotOwn(data1, data2, par)
ind_correctDir = [3 2 1 8 7 6 5 4]; 

if contains(par.InField, {'AE', 'angular','Hasson'}, 'IgnoreCase', true)
    data1_rad = abs(data1(1,ind_correctDir)); 
    data2_rad = abs(data2(1,ind_correctDir));
else
    data1_rad = data1(1,ind_correctDir);
    data2_rad = data2(1,ind_correctDir);   
end

% Gestione NaN
if strcmp(par.mod, 'Zero')
    data1_rad(isnan(data1_rad)) = 0;
    data2_rad(isnan(data2_rad)) = 0;
elseif strcmp(par.mod, 'Interp')
    if any(isnan(data1_rad))
        data1_rad = fillmissing(data1_rad, 'linear', 'EndValues', 'nearest');
    end
    if any(isnan(data2_rad))
        data2_rad = fillmissing(data2_rad, 'linear', 'EndValues', 'nearest');
    end
end

max_v = par.maxValue;
if isfield(par, 'minValue')
    min_v = par.minValue;
else
    min_v = 0;
end

if min_v< 0
    min_v = 0;
end

theta = linspace(0, 2*pi, 9);
% Prepara poligoni chiusi
p1 = [data1_rad, data1_rad(1)];
p2 = [data2_rad, data2_rad(1)];

figure('Name', par.figureName, 'Color', 'w', 'Position', [100 100 800 700]);
nexttile;

% === FONT FLESSIBILI (IDENTICO AGLI ALTRI) ===
fs_ticks_default = 10;
fs_legend_default = 10;
if isfield(par, 'font_ticks'),  fs_ticks = par.font_ticks; else, fs_ticks = fs_ticks_default; end
if isfield(par, 'font_legend'), fs_legend = par.font_legend; else, fs_legend = fs_legend_default; end
if isfield(par, 'font_title'),  fs_title_mult = par.font_title / fs_ticks; else, fs_title_mult = 1.2; end
% === FINE FONT ===

lineStyle1 = par.lineStyle1; color1 = par.color1;
h1 = polarplot(theta, p1, 'LineStyle', lineStyle1, 'LineWidth', 2, ...
    'Color', color1, 'HandleVisibility', 'on'); hold on;
polarplot(theta, p1, 'o', 'MarkerFaceColor', color1, ...
    'MarkerEdgeColor', color1, 'HandleVisibility', 'off', 'MarkerSize', 4);

lineStyle2 = par.lineStyle2; color2 = par.color2;
h2 = polarplot(theta, p2, 'LineStyle', lineStyle2, 'LineWidth', 2, ...
    'Color', color2, 'HandleVisibility', 'on'); hold on;
polarplot(theta, p2, 'o', 'MarkerFaceColor', color2, ...
    'MarkerEdgeColor', color2, 'HandleVisibility', 'off', 'MarkerSize', 4);

% Stelle p-value
if isfield(par, 'p_vals') && ~isempty(par.p_vals)
    p_vals_curr = par.p_vals(ind_correctDir);
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

% Titoli corretti
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
% Assi + FONT (sostituito solo FontSize)
pax.ThetaTick = linspace(0, 315, 8);
pax.ThetaTickLabel = {'D3','D2','D1','D8','D7','D6','D5','D4'};
pax.FontWeight = 'bold'; 
pax.FontSize = fs_ticks;                    % ← Override
pax.TitleFontSizeMultiplier = fs_title_mult; % ← Override
pax.RLim = [min_v max_v];                   % ← Zoom

% Leggenda con override
legendNames = par.legendNames;
legendNames_plot = cellfun(@(c) strjoin(c, '\n '), legendNames, 'UniformOutput', false);
lgd = legend([h1 h2], legendNames_plot, 'Orientation', 'horizontal', 'Interpreter', 'tex');
lgd.FontSize = fs_legend;  % ← Override
lgd.Layout.Tile = 'south';

set(gcf, 'WindowState', 'maximized');
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])

% Salvataggio (invariato)
pMaxDistDiff_path = [par.save_dir,filesep];
par.savePlotEpsPdfMat.dir_png = fullfile(pMaxDistDiff_path, 'PNGs');
par.savePlotEpsPdfMat.dir_pdf = fullfile(pMaxDistDiff_path, 'PDFs');
par.savePlotEpsPdfMat.dir_mat = fullfile(pMaxDistDiff_path, 'MATfiles');
name_fig =['PolarPlot_',par.InField,'_',par.figureName];
par.savePlotEpsPdfMat.file_name = name_fig;
savePlotEpsPdfMat(gcf,par.savePlotEpsPdfMat)
close(gcf)
end
