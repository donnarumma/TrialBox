function outMatrix = structToMatrix(singleStruct, fieldName)
% structToMatrix - Converte un campo della struct in una matrice 3x8
%
% INPUT:
%   singleStruct: la tua stats.single (1x24)
%   fieldName: stringa col nome del campo (es. 'maxDist', 'RT', 'timeMax')
%
% OUTPUT:
%   outMatrix: matrice 3x8 (righe=condizioni, colonne=direzioni)

% Inizializziamo la matrice con NaN per sicurezza
outMatrix = nan(3, 8);

for i = 1:numel(singleStruct)
    % Recuperiamo riga (cond) e colonna (dir)
    try
        r = singleStruct(i).cond;
        c = singleStruct(i).dir;
    catch
        r = singleStruct(i).Condition;
        c = singleStruct(i).Direction;
    end
    
    % Assegniamo il valore del campo richiesto alla cella corrispondente
    % Usiamo il dynamic field indexing singleStruct(i).(fieldName)
    outMatrix(r, c) = singleStruct(i).(fieldName);
end

end