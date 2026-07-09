function savePlotEpsPdfMat(gcf,par)


if ~exist(par.dir_png, 'dir')
    mkdir(par.dir_png);
end

if ~exist(par.dir_pdf, 'dir')
    mkdir(par.dir_pdf);
end

if ~exist(par.dir_mat, 'dir')
    mkdir(par.dir_mat);
end

file_name_png = strcat(par.file_name, '.png');
file_name_pdf = strcat(par.file_name, '.pdf');
file_name_mat = strcat(par.file_name, '.mat');


% PNG
set(gcf, 'InvertHardcopy', 'off');
fulldir_png = fullfile(par.dir_png, file_name_png);
saveas(gcf, fulldir_png, 'png');

% PDF
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'PaperPositionMode', 'auto'); % Usa le dimensioni della figura per il PDF
% set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Assicurati che la figura occupi tutto lo spazio disponibile
set(gcf, 'PaperOrientation', 'landscape');

fulldir_pdf = fullfile(par.dir_pdf, file_name_pdf);
% saveas(gcf, fulldir_pdf, 'pdf');
exportgraphics(gcf, fulldir_pdf, 'ContentType', 'vector'); % Salva come PDF vettoriale

% MAT
fulldir_fig = fullfile(par.dir_mat, strrep(file_name_mat, '.mat', '.fig'));
savefig(gcf, fulldir_fig);


% % EPS
% set(gcf, 'InvertHardcopy', 'off');
% fulldir_eps = fullfile(par.dir_eps, file_name_eps);
% saveas(gcf, fulldir_eps, 'epsc');
% % print(gcf, fulldir_eps, '-depsc', '-r300');
