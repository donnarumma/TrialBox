function Var = extractVarParams(S,par)
    % S: la tua struct (RT_trial, PV_trial, PVT_trial)
    % paramType: 'ms' per RT/PVT, 'NaN' for PV
    
    num_dir = 8;
    Var = struct();

    for d = 1:num_dir
        dir_name = sprintf('dir%d', d);
        
        % Estrazione dei vettori 24x1
        soloS  = S(1).(dir_name); % Riga 1
        soloK  = S(2).(dir_name); % Riga 2
        jointS = S(3).(dir_name); % Riga 3
        jointK = S(4).(dir_name); % Riga 4
        
        if strcmp(par.time,'ms')
            Var.S.Solo(:,d)     = 1000*mean(soloS);
            Var.S.Joint(:,d)    = 1000*mean(jointS);
            Var.K.Solo(:,d)     = 1000*mean(soloK);
            Var.K.Joint(:,d)    = 1000*mean(jointK);
        else
            Var.S.Solo(:,d)     = mean(soloS);
            Var.S.Joint(:,d)    = mean(jointS);
            Var.K.Solo(:,d)     = mean(soloK);
            Var.K.Joint(:,d)    = mean(jointK);
        end
    end
end