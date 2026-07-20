function [DeltaJoint, DeltaSolo,Var] = extractBehavDeltas(Jxy_S, Jxy_K, paramName)
    num_sess = length(Jxy_K);
    num_dir = 8;
    DeltaJoint = zeros(num_sess, num_dir);
    DeltaSolo  = zeros(num_sess, num_dir);
    Var = struct();
    for s = 1:num_sess
        for d = 1:num_dir
            dir_name = sprintf('dir%d', d);
            
            % --- ESTRAZIONE DATI ---
            % Assumiamo: riga 1=SoloS, riga 2=SoloK, riga 3=Joint
            valS_Solo  = mean(Jxy_S(s).(paramName)(1).(dir_name), 'omitnan');
            valK_Solo  = mean(Jxy_K(s).(paramName)(2).(dir_name), 'omitnan');
            
            valS_Joint = mean(Jxy_S(s).(paramName)(3).(dir_name), 'omitnan');
            valK_Joint = mean(Jxy_K(s).(paramName)(3).(dir_name), 'omitnan');

            % --- CALCOLO DELTA (K - S per parametri "meno è meglio") ---
            % Delta Solo: Rappresenta chi è più bravo "di natura"
            DeltaSolo(s,d) = valK_Solo - valS_Solo;
            
            % Delta Joint: Rappresenta chi è più bravo "in coppia"
            DeltaJoint(s,d) = valK_Joint - valS_Joint;
            Var.S.Solo(s,d)     = valS_Solo;
            Var.S.Joint(s,d)    = valS_Joint;
            Var.K.Solo(s,d)     = valK_Solo;
            Var.K.Joint(s,d)    = valK_Joint;
        end
    end
end