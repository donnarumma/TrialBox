function confError = confusionError(cmatrix)

% Dimensione della matrice cmatrix
[n_cmatrix, ~] = size(cmatrix);

% Inizializzamo la struttura confError
confError = struct();

% Inizializzamo l'indice della struttura
index = 1;
% Ciclo per creare gli elementi della struttura confError
for i = 1:n_cmatrix
    for j = i+1:n_cmatrix % Scorre solo la metà superiore (escludendo la diagonale principale)
        confError(index).class_i = i;
        confError(index).class_j = j;
        confError(index).sum_ij = cmatrix(i, j) + cmatrix(j, i);
        index = index + 1;
    end
end

[~,indMax] = max([confError.sum_ij]);
disp(confError(indMax));