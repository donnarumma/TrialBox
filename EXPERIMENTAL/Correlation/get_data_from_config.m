function [active_neural, active_label, active_embed, ...
          obs_neural, obs_label, obs_embed] = get_data_from_config(config, DATA,new_indices)

cond_key = sprintf('cond_%d', config.condition);
subject_a = config.active_subject;
subject_o = config.obs_subject;

% Se new_indices non è fornito o è vuoto, prendi tutto
if nargin < 3 || isempty(new_indices)
    active_idx = ':';
    obs_idx    = ':';
else
    active_idx = new_indices;
    obs_idx    = new_indices;
end

if config.condition < 3
    % Condizione 1 o 2: uso i campi con suffisso _active e _obs
    active_neural  = DATA.(subject_a).(cond_key).neural_active(active_idx , :);
    active_label   = DATA.(subject_a).(cond_key).label_active(active_idx );
    active_embed   = DATA.(subject_a).(cond_key).embed_active(active_idx , :);

    obs_neural = DATA.(subject_o).(cond_key).neural_obs(obs_idx, :);
    obs_label  = DATA.(subject_o).(cond_key).label_obs(obs_idx);
    obs_embed  = DATA.(subject_o).(cond_key).embed_obs(obs_idx, :);

elseif config.condition == 3 
    % Condizione 3 (dati joint): uso nomi neutri (neural, label, embed)
    active_neural  = DATA.(subject_a).(cond_key).neural(active_idx, :);
    active_label   = DATA.(subject_a).(cond_key).label(active_idx);
    active_embed   = DATA.(subject_a).(cond_key).embed(active_idx, :);

    obs_neural = DATA.(subject_o).(cond_key).neural(obs_idx, :);
    obs_label  = DATA.(subject_o).(cond_key).label(obs_idx);
    obs_embed  = DATA.(subject_o).(cond_key).embed(obs_idx, :);
%elseif config.condition == 4 
   % continue
    
end