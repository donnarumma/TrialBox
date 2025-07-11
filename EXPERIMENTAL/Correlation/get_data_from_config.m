function [active_neural, active_label, active_embed, ...
          obs_neural, obs_label, obs_embed] = get_data_from_config(config, DATA, new_indices)

cond_key = sprintf('cond_%d', config.condition);
subject_a = config.active_subject;
subject_o = config.obs_subject;

if config.condition < 3
    % Condizione 1 o 2: uso i campi con suffisso _active e _obs
    active_neural  = DATA.(subject_a).(cond_key).neural_active(new_indices, :);
    active_label   = DATA.(subject_a).(cond_key).label_active(new_indices);
    active_embed   = DATA.(subject_a).(cond_key).embed_active(new_indices, :);

    obs_neural = DATA.(subject_o).(cond_key).neural_obs(new_indices, :);
    obs_label  = DATA.(subject_o).(cond_key).label_obs(new_indices);
    obs_embed  = DATA.(subject_o).(cond_key).embed_obs(new_indices, :);

elseif config.condition == 3 
    % Condizione 3 (dati joint): uso nomi neutri (neural, label, embed)
    active_neural  = DATA.(subject_a).(cond_key).neural(new_indices, :);
    active_label   = DATA.(subject_a).(cond_key).label(new_indices);
    active_embed   = DATA.(subject_a).(cond_key).embed(new_indices, :);

    obs_neural = DATA.(subject_o).(cond_key).neural(new_indices, :);
    obs_label  = DATA.(subject_o).(cond_key).label(new_indices);
    obs_embed  = DATA.(subject_o).(cond_key).embed(new_indices, :);
elseif config.condition == 4 
    continue
    

end