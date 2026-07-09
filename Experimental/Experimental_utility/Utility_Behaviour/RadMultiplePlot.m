function RadMultiplePlot(varargin)
% RADMULTIPLEPLOT - Polarplot multi-curva (1-4 linee) custom
% RadMultiplePlot(data1,data2,data3,data4,par)
% Es: RadMultiplePlot(MNDT_ActK,MNDT_ObsS,MNDT_JointK,par.RadMultiplePlot)
data_cell = varargin(1:end-1);
par = varargin{end};
ind_correctDir = [3 2 1 8 7 6 5 4];

% Preallocazione fissa MAX 4 curve
data_rad = cell(1,4);
n_curves = 0;
% Filtra dati validi
for i = 1:length(data_cell)
    if ~isempty(data_cell{i}) && ~all(isnan(data_cell{i}))
        n_curves = n_curves + 1;
        data_rad{n_curves} = data_cell{i}(ind_correctDir);
    end
end
if n_curves == 0
    warning('Nessun dato valido!'); return;
end
if isfield(par, 'InField') && contains(par.InField, {'AE', 'angular', 'Hasson'}, 'IgnoreCase', true)
    for i = 1:n_curves
        data_rad{i} = abs(data_rad{i});
    end
end
% Gestione NaN
if isfield(par, 'mod')
    for i = 1:n_curves
        if strcmp(par.mod, 'Zero')
            data_rad{i}(isnan(data_rad{i})) = 0;
        elseif strcmp(par.mod, 'Interp')
            data_rad{i} = fillmissing(data_rad{i}, 'linear', 'EndValues', 'nearest');
        end
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

figure('Name', par.figureName, 'Color', 'w', 'Position', [100 100 1000 800]);
ax = polaraxes;
hold on;

% === FONT FLESSIBILI (DEFAULT + OVERRIDE) ===
fs_ticks_default = 10;
fs_legend_default = 10;
if isfield(par, 'font_ticks'),  fs_ticks = par.font_ticks; else, fs_ticks = fs_ticks_default; end
if isfield(par, 'font_legend'), fs_legend = par.font_legend; else, fs_legend = fs_legend_default; end
if isfield(par, 'font_title'),  fs_title_mult = par.font_title / fs_ticks; else, fs_title_mult = 1.2; end

ax.FontSize = fs_ticks;
ax.TitleFontSizeMultiplier = fs_title_mult;
ax.FontWeight = 'bold';
% === FINE FONT ===

handles = gobjects(1, n_curves);

% === COLORS E STYLES CORRETTI (sicuri) ===
colors = {par.color1, par.color2, par.color3, par.color4};
styles = {par.lineStyle1, par.lineStyle2, par.lineStyle3, par.lineStyle4};

% Plot curve
for i = 1:n_curves
    p = [data_rad{i}, data_rad{i}(1)];  % Poligono chiuso
    
    handles(i) = polarplot(ax, theta, p, 'LineStyle', styles{i}, ...
        'LineWidth', 2, 'Color', colors{i}, 'HandleVisibility', 'on'); 
    hold(ax, 'on');
    
    % Marker
    polarplot(ax, theta(1:8), data_rad{i}, 'o', ...
        'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', colors{i}, ...
        'HandleVisibility', 'off', 'MarkerSize', 4);
end

% Stelle p-value (font fisso 18 OK)
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

% Titoli colorati
titles = {};
if isfield(par, 'cond_names1') && ~isempty(par.cond_names1) && ~isempty(par.cond_names1{1})
    colorStr1 = sprintf('%.2g,%.2g,%.2g', par.color1);
    titles{end+1} = ['\color[rgb]{', colorStr1, '}', par.cond_names1{1}];
end
if isfield(par, 'cond_names2') && ~isempty(par.cond_names2) && ~isempty(par.cond_names2{1})
    colorStr2 = sprintf('%.2g,%.2g,%.2g', par.color2);
    titles{end+1} = ['\color[rgb]{', colorStr2, '}', par.cond_names2{1}];
end
if isfield(par, 'cond_names3') && ~isempty(par.cond_names3) && ~isempty(par.cond_names3{1})
    colorStr3 = sprintf('%.2g,%.2g,%.2g', par.color3);
    titles{end+1} = ['\color[rgb]{', colorStr3, '}', par.cond_names3{1}];
end
if isfield(par, 'cond_names4') && ~isempty(par.cond_names4) && ~isempty(par.cond_names4{1})
    colorStr4 = sprintf('%.2g,%.2g,%.2g', par.color4);
    titles{end+1} = ['\color[rgb]{', colorStr4, '}', par.cond_names4{1}];
end
title(ax, titles, 'Interpreter', 'tex');

% Assi polari
ax.ThetaTick = linspace(0, 315, 8);
ax.ThetaTickLabel = {'D3','D2','D1','D8','D7','D6','D5','D4'};
if isfield(par, 'unitMeasure')
    ax.RAxis.Label.String = par.unitMeasure;
    ax.RAxis.Label.Interpreter = 'tex';
end
ax.RLim = [min_v max_v];

% Legenda dinamica
valid_legends = par.legendNames(1:n_curves);
legendNames_plot = cellfun(@(c) strjoin(c, '\n '), valid_legends, 'UniformOutput', false);
lgd = legend(ax, handles, legendNames_plot, 'Orientation', 'horizontal', 'Location', 'southoutside');
lgd.Interpreter = 'tex';
lgd.FontSize = fs_legend;  % Usa default o override!

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
