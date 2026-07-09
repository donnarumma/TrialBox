function [Jxy_S, Jxy_K] = JxyMatrix(Result)
    monkey_names = {'S', 'K'};
    fields_to_extract = {'exit_time', 'angular_error', 'is_correct','CMT'};
    data_out = struct();

    for sj = 1:length(monkey_names)
        m_name = monkey_names{sj};
        struct_name = ['Jxy_', m_name];
        
        for s = 1:length(Result)
            data_out.(struct_name)(s).session = Result(s).session;
            data_out.(struct_name)(s).monkey  = m_name;
            
            for cond = 1:3
                for d = 1:8
                    dirName = sprintf('dir%d', d);
                    trials = Result(s).(m_name)(cond).(dirName);
                    
                    if ~isempty(trials)
                        num_trials = length(trials);
                        temp_data = struct();
                        for f = 1:length(fields_to_extract)
                            temp_data.(fields_to_extract{f}) = NaN(1, num_trials);
                        end
                        
                        for t = 1:num_trials
                            curr_info = trials(t).info;
                            for f = 1:length(fields_to_extract)
                                fname = fields_to_extract{f};
                                if isfield(curr_info, fname) && ~isempty(curr_info.(fname))
                                    temp_data.(fname)(t) = curr_info.(fname);
                                end
                            end
                        end
                        
                        % --- ASSEGNAZIONE VETTORI ---
                        for f = 1:length(fields_to_extract)
                            fname = fields_to_extract{f};
                            data_out.(struct_name)(s).(fname)(cond).(dirName) = temp_data.(fname);
                        end
                        
                        % --- AGGIUNTA EXIT PERCENTAGE (Scalare) ---
                        % Calcoliamo la percentuale di 1 rispetto ai trial validi
                        % valid_trials = temp_data.is_correct(~isnan(temp_data.is_correct));
                        % if ~isempty(valid_trials)
                        %     data_out.(struct_name)(s).ExitCorr(cond).(dirName) = mean(valid_trials) * 100;
                        % else
                        %     data_out.(struct_name)(s).ExitCorr(cond).(dirName) = 0;
                        % end
                        if any(isnan(temp_data.exit_time))
                            % Se c'è almeno un NaN, la percentuale diventa NaN
                            data_out.(struct_name)(s).ExitCorr(cond).(dirName) = NaN;
                            data_out.(struct_name)(s).CMT(cond).(dirName) = NaN;                            
                        else
                            % Altrimenti (tutti i trial sono validi), calcoliamo la media
                            % Usiamo comunque is_correct per determinare il successo
                            data_out.(struct_name)(s).ExitCorr(cond).(dirName) = mean(temp_data.is_correct) * 100;
                        end
                    else
                        % Caso direzione vuota
                        for f = 1:length(fields_to_extract)
                            data_out.(struct_name)(s).(fields_to_extract{f})(cond).(dirName) = [];
                        end
                        data_out.(struct_name)(s).ExitCorr(cond).(dirName) = NaN;
                    end
                end
            end
        end
    end
    
    Jxy_S = data_out.Jxy_S;
    Jxy_K = data_out.Jxy_K;
end