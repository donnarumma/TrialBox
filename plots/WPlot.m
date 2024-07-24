function [W_out,MI_Filter_ID] = WPlot(sorted_IndMI,Wcsp)

% Calcolo del numero di righe e colonne
num_rows = size(Wcsp{1, 1},2);
num_cols = size(Wcsp,2);
n_class = size(Wcsp,1);
nch = size(Wcsp{1, 1},1);

% Numero totale di elementi
n_features = num_rows * num_cols;

% Creare un vettore con numeri in ordine crescente da 1 a n_features
elements = 1:n_features;

% Reshape il vettore in una matrice con num_rows righe e num_cols colonne
Filter_Matrix = reshape(elements, num_rows, num_cols);

% Inizializza una struttura per memorizzare gli indici delle colonne e le posizioni
MI_Filter_ID = struct('Filter', cell(size(sorted_IndMI)), 'PosFilter', cell(size(sorted_IndMI)));

% Trova la colonna corrispondente e la posizione per ciascun valore in sorted_IndMI
for i = 1:numel(sorted_IndMI)
    values = sorted_IndMI{i};
    column_indices = zeros(size(values));
    position_indices = zeros(size(values));
    
    for j = 1:numel(values)
        [row, col] = find(Filter_Matrix == values(j));
         column_indices(j) = col;
         position_indices(j) = row;
    end
    
    MI_Filter_ID(i).Filter = column_indices;
    MI_Filter_ID(i).PosFilter = position_indices;
end


% Inizializza una cella per memorizzare l'output
W_out = cell(n_class, 1);

% Itera attraverso ciascun elemento di Filter_ID
for i = 1:n_class
    filter_indices = MI_Filter_ID(i).Filter;
    pos_indices = MI_Filter_ID(i).PosFilter;
    
    % Prealloca la matrice per i dati estratti
    extracted_data = zeros(nch, numel(filter_indices));
    
    % Itera attraverso ciascun filtro e posizione
    for j = 1:numel(filter_indices)
        filter_idx = filter_indices(j);
        pos_idx = pos_indices(j);
        
        % Estrai i dati dalla matrice corrispondente in W
        data_app = Wcsp{i, filter_idx}(:, pos_idx);
        
        % Inserisci i dati estratti nella posizione corretta
        extracted_data(:, j) = data_app;
    end
    
    % Memorizza i dati estratti in W_out
    W_out{i} = extracted_data;
end
