function TopoplotMi(W1,W2,par)

addpath(genpath("D:\eeglab2023.1\"))
load("D:\main_scriptNSA\selected_chanlocs.mat");

conditionTitles = {'Left Hand', 'Right Hand', 'Foot', 'Tongue'};


% Creazione di una figura
figure
for i = 1:size(W1,1)
    % Crea una subplot per la condizione corrente
    subplot(2, 4, i);
    W1plot = W1{i,1};
    % Plot topografico utilizzando 'topoplot'
    topoplot(W1plot(:,1), selected_chanlocs, 'maplimits', 'absmax', 'electrodes','pts','style', 'both', 'emarkersize', 10);
    % 'ptslabels');
    colorbar;
    title(sprintf( 'Accuracy: %.1f %.s\n', par.accuracy1(i)),conditionTitles{i}, 'FontSize', 14);

    hold on
    subplot(2, 4, i+4);
    W2plot = W2{i,1};
    % Plot topografico utilizzando 'topoplot'
    topoplot(W2plot(:,1), selected_chanlocs, 'maplimits', 'absmax', 'electrodes','pts','style', 'both', 'emarkersize', 10);
    % 'ptslabels');
    colorbar;
    title(sprintf( 'Accuracy: %.1f %.s\n', par.accuracy2(i)),conditionTitles{i}, 'FontSize', 14);
end

sgtitle(par.titleAll);

% hColorbar = colorbar('southoutside');
% pos = get(hColorbar, 'Position');
% set(hColorbar, 'Position', [pos(1), pos(2)-0.05, pos(3), pos(4)]); % Aumenta lo spazio tra i subplot e la colorbar
% Aggiungere titoli al centro di ogni riga di subplot
% Calcolare le posizioni per i titoli centrali
xPos = 0.5; % Posizione centrale lungo l'asse x
yPos1 = 0.95; % Posizione per il titolo della prima riga (sopra i subplot)
yPos2 = 0.45; % Posizione per il titolo della seconda riga (sopra i subplot della seconda riga)

% Aggiungere i titoli con la funzione text
annotation('textbox', [xPos - 0.25, yPos1, 0.5, 0.05], 'String', par.title1, 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 14, 'FontWeight', 'bold');
annotation('textbox', [xPos - 0.25, yPos2, 0.5, 0.05], 'String', par.title2, 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 14, 'FontWeight', 'bold');
