function outMatrix = meanTrialsToMatrix(trialStruct, fieldName, num_cond, num_dir)
% MEANTRIALSTOMATRIX - Medie per cond×dir DA struct di SINGOLI TRIAL
%
% INPUT:
%   trialStruct: 1xN_trials struct da compute_pca_distance_params_trials().trial
%   fieldName:   'maxDist', 'RT', ecc.
%   num_cond:    3 (default)
%   num_dir:     8 (default)
%
% OUTPUT:
%   outMatrix: [num_cond × num_dir] medie dei valori fieldName

if nargin < 3, num_cond = 3; end
if nargin < 4, num_dir  = 8; end

outMatrix = nan(num_cond, num_dir);

% Accumulatore per conteggi (debug)
trialCount = zeros(num_cond, num_dir);

for i = 1:numel(trialStruct)
    r = trialStruct(i).cond;
    c = trialStruct(i).dir;
    
    % Somma incrementale (per media)
    if isnan(outMatrix(r,c))
        outMatrix(r,c) = trialStruct(i).(fieldName);
    else
        outMatrix(r,c) = outMatrix(r,c) + trialStruct(i).(fieldName);
    end
    trialCount(r,c) = trialCount(r,c) + 1;
end

% Divisione per numero trial per cella
for r = 1:num_cond
    for c = 1:num_dir
        if trialCount(r,c) > 0
            outMatrix(r,c) = outMatrix(r,c) / trialCount(r,c);
        end
    end
end
end
