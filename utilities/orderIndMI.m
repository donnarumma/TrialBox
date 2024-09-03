function sorted_IndMI = orderIndMI(IndMI,Pos_IndMI)

% Numero di elementi nelle celle
num_cells = numel(IndMI);

% Inizializza una cella per i valori ordinati
sorted_IndMI = cell(num_cells, 1);

% Ciclo per ordinare i valori in base agli indici
for i = 1:num_cells
    values = IndMI{i};
    indices = Pos_IndMI{i};
    
    % Verifica che il numero di valori e indici sia uguale
    if numel(values) == numel(indices)
        % Crea un vettore vuoto per i valori ordinati
        sorted_values = zeros(size(values));
        
        % Popola il vettore ordinato in base agli indici
        sorted_values(indices) = values;
        
        % Assegna i valori ordinati alla cella di output
        sorted_IndMI{i} = sorted_values;
    else
        error('Il numero di valori e indici non corrisponde nella cella %d', i);
    end
end