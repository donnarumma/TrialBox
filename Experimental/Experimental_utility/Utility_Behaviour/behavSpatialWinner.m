%function Table = behavSpatialWinner(data, par)

function Table = behavSpatialWinner(data, par)
    % BEHAVSPATIALWINNER Valuta i vincitori per Bias, Shift e Reliability
    % basandosi su una soglia fisica di significatività (gradi).
    
    num_dir = par.num_dir;
    
    % --- CONFIGURAZIONE SOGLIA FISICA ---
    % Se il settore è 45°, 2.25° è il 5% dello spazio totale.
    % Differenze sotto questa soglia sono considerate "Balanced".
    deg_threshold = 2.25; 
    rel_percent_threshold = 0.10; % Soglia 10% per la Reliability (Relativa)
    % Definizione struttura tabella
    VarNames = {'Direction','Bias_Solo','Bias_Joint','Bias_S','Bias_K','Shift','Reliability'};
    VarTypes = {'double','string','string','string','string','string','string'};
    Table = table('Size',[num_dir, length(VarNames)], 'VariableNames', VarNames, 'VariableTypes', VarTypes);

    % Metadati per la visualizzazione
    Table.Properties.VariableUnits = {'', 'Chi si allontana meno', 'Chi si allontana meno', ...
        'Riduce il bias', 'Riduce il bias', 'Non varia il comportamento', 'è più costante'};

    for ndir = 1:num_dir
        % Estrazione dati Bias
        BS_Slo = data.BiasS_Solo(ndir);
        BS_Jnt = data.BiasS_Joint(ndir);
        BK_Slo = data.BiasK_Solo(ndir);
        BK_Jnt = data.BiasK_Joint(ndir);

        % --- LOGICA BIAS (Puntuale) ---
        % Chi ha l'errore sistematico minore?
        if BS_Slo < BK_Slo, WB_solo = 'S'; else, WB_solo = 'K'; end
        if BS_Jnt < BK_Jnt, WB_joint = 'S'; else, WB_joint = 'K'; end
        
        % Il Joint ha migliorato (Correction) o peggiorato (Influence) il singolo?
        if BS_Jnt < BS_Slo, wB_s = 'Correction'; else, wB_s = 'Influence'; end
        if BK_Jnt < BK_Slo, wB_k = 'Correction'; else, wB_k = 'Influence'; end

        % --- LOGICA SHIFT (Soglia Fisica) ---
        % Chi è il Leader (chi ha variato meno la traiettoria tra Solo e Joint)?
        ShiftS = data.Shift_S(ndir);
        ShiftK = data.Shift_K(ndir);
        diffSh = ShiftS - ShiftK; 

        if abs(diffSh) < deg_threshold
            sh = 'Balanced';
        elseif diffSh < 0  % S ha cambiato meno gradi di K -> S è Leader
            sh = 'S';
        else               % K ha cambiato meno gradi di S -> K è Leader
            sh = 'K';
        end

        % --- LOGICA RELIABILITY (Soglia Fisica) ---
        % Chi è più costante (ha meno dispersione in gradi)?
        % Nota: assumendo che Reliability qui sia espressa come variabilità (gradi)
        RelS = data.Reliability_S(ndir);
        RelK = data.Reliability_K(ndir);
        diffRel_abs = abs(RelS - RelK);
        maxRel = max(abs(RelS), abs(RelK));
        
        % Se la differenza è meno del 10% del valore massimo, sono Balanced
        if diffRel_abs < (maxRel * rel_percent_threshold)
            Rel = 'Balanced';
        elseif RelS < RelK  % Assumendo meno dispersione = migliore
            Rel = 'S';
        else
            Rel = 'K';
        end

        % Riempimento riga tabella
        Table(ndir,:) = {ndir, WB_solo, WB_joint, wB_s, wB_k, sh, Rel};
    end
end