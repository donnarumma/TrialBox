function [Ind_sel,position_Ind,ind_com,Entr_sel] = MutualInformationComplete(data,m,k)

% function Ind_final = MIBIF_decode(V,label,m,k)
% MIBIV: Mutual Information-Based Best Individual Feature (MIBIV)
% This function evaluates the index of the best features for each class.
% INPUT:
% 'data': a matrix containing on the row the trials and on the column the
% features
% 'label': a vector containing the class labels
% 'm': the number elements to consider
% 'k': the maximum number of mutual information elements to include as features
% OUTPUT:
% 'Ind_final': the final index containing the indices of the projection matrix features to consider for classification.


numTrial = size(data,1);
len_Vi = size(data,2); % Numero di features

% Entropia globale (senza separazione per classe)
P = mean(label); % Probabilità media delle classi
H = -sum(P .* log2(P + eps)); % Entropia totale
I = zeros(len_Vi,1); % MI per ogni feature

for j = 1:len_Vi
    Hc = 0; % Entropia condizionale

    for i = 1:numTrial
        % Parzen Window
        sig = std(data(i,j) - data(:,j));
        h = sig * (4 / (3 * numTrial))^(1/5);
        phi = exp(-((data(i,j) - data(:,j)).^2) / (2 * h^2)) / sqrt(2 * pi);

        % Stima della probabilità condizionale senza usare `sel_class`
        cp = sum(phi(label == 1)) / sum(label == 1);
        np = sum(phi(label == 0)) / sum(label == 0);

        % Probabilità della feature j
        p = cp * P + np * (1 - P);

        % Probabilità condizionali
        pc = cp * P / (p + eps);
        pn = np * (1 - P) / (p + eps);

        % Entropia condizionale
        Hc = Hc + (pc * log2(pc + eps) + pn * log2(pn + eps));
    end

    % Mutual Information della feature j
    I(j) = H + Hc; % Ora indipendente dalla classe
end

% Selezione delle feature
[Entr, ind] = sort(I, 'descend'); % Ordine decrescente
ind_sel = ind(1:k); % Prime k features migliori

% Indici complementari (CSP)
if (m == 0)
    ind_com = ind_sel;
else
    ind_com = (4*m) * ceil(ind_sel / (2*m)) - (2*m - 1) - ind_sel;
end

% Indici finali
[Ind_sel, position_Ind] = unique([ind_sel ind_com]);
Entr_sel = Entr(Ind_sel);

