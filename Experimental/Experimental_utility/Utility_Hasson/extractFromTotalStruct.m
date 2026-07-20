function [DeltaJoint, DeltaSolo,Var] = extractFromTotalStruct(S, paramType)
    % S: la tua struct (RT_total, VM_total, ecc.)
    % paramType: 'less' per RT/T(VM), 'more' for VM
    
    num_sess = length(S(1).dir1);
    num_dir = 8;
    DeltaJoint = zeros(num_sess, num_dir);
    DeltaSolo  = zeros(num_sess, num_dir);
    Var = struct();

    for d = 1:num_dir
        dir_name = sprintf('dir%d', d);
        
        % Estrazione dei vettori 24x1
        soloS  = S(1).(dir_name); % Riga 1
        soloK  = S(2).(dir_name); % Riga 2
        jointS = S(3).(dir_name); % Riga 3
        jointK = S(4).(dir_name); % Riga 4

        if strcmp(paramType, 'less')
            % Logica K - S (RT, T(VM))
            DeltaSolo(:, d)  = soloK - soloS;
            DeltaJoint(:, d) = jointK - jointS;
        else
            % Logica S - K (VM)
            DeltaSolo(:, d)  = soloS - soloK;
            DeltaJoint(:, d) = jointS - jointK;
        end
        Var.S.Solo(:,d)     = soloS;
        Var.S.Joint(:,d)    = jointS;
        Var.K.Solo(:,d)     = soloK;
        Var.K.Joint(:,d)    = jointK;
    end
end