function hfig = SubGridPlot(data,par)


% Creazione della griglia di plot
hfig = figure;
imagesc(data); % Visualizza la matrice con colori
colormap(jet); % Utilizza la scala di colori jet (puoi usare 'gray' per la scala di grigi)
% colormap("gray");
colorbar; % Aggiunge una barra dei colori per riferimento
axis square; % Mantiene l'asse quadrato
% Aggiungi il titolo e spostalo più in alto possibile
t = title(par.title);
set(t, 'Units', 'normalized'); % Imposta l'unità a 'normalized' per un posizionamento più facile
titlePosition = get(t, 'Position');
titlePosition(2) = titlePosition(2) + 0.05; % Sposta il titolo verso l'alto
set(t, 'Position', titlePosition);
% Aggiunge i valori all'interno delle celle
for i = 1:9
    for j = 1:9
        text(j, i, num2str(data(i, j), '%.2f'), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
end

% Aggiunge le coordinate al bordo della griglia
set(gca, 'XTick', 1:9, 'XTickLabel', 1:9, 'YTick', 1:9, 'YTickLabel', 1:9);
axis off
% Aggiunge le etichette delle combinazioni (facoltativo)
for i = 1:9
    for j = 1:9
        % Aggiungi le coordinate sopra e a sinistra della griglia
        text(j, 0, num2str(j), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        text(0.1, i, num2str(i), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
end

% t = title(par.title);
% set(t, 'Position', get(t, 'Position') + [0, 1, 0]); % Sposta il titolo verso l'alto
% xlabel('Test');
% ylabel('Train');
