function Data = build24BehaviourICDStruct(A)
% BUILD24BEHAVIOURICDSTRUCT Convert ICD matrices into a 24-entry structure.
%
% This function converts the Inter-Cursor Distance (ICD) metrics stored in
% condition-by-direction matrices into a standardized 1x24 behavioral
% structure, matching the format used throughout the behavioral pipeline.
%
% INPUT
% -------------------------------------------------------------------------
% A : struct
%
%   Structure containing ICD summary metrics:
%
%       A.meanMatrix   [3 x 8]
%           Mean ICD values.
%
%       A.maxMatrix    [3 x 8]
%           Maximum ICD values.
%
%       A.aucMatrix    [3 x 8]
%           Area-under-the-curve (AUC) ICD values.
%
%   Rows correspond to experimental conditions:
%       1 = Solo S
%       2 = Solo K
%       3 = Joint
%
%   Columns correspond to movement directions (1-8).
%
% OUTPUT
% -------------------------------------------------------------------------
% Data : struct (1 x 24)
%
%   Each element corresponds to one condition-direction pair:
%
%       Data(i).ICD_mean
%       Data(i).ICD_max
%       Data(i).ICD_auc
%
%   The ordering follows the convention adopted in the behavioral
%   analyses:
%
%       [Cond1 Dir1-8, Cond2 Dir1-8, Cond3 Dir1-8]
%
%   For example:
%
%       Data(1:8)    -> Condition 1
%       Data(9:16)   -> Condition 2
%       Data(17:24)  -> Condition 3
%
% Author: Mirco Frosolone
% -------------------------------------------------------------------------

%% Validate input fields

requiredFields = {'meanMatrix','maxMatrix','aucMatrix'};

for iField = 1:numel(requiredFields)
    if ~isfield(A, requiredFields{iField})
        error('Missing required field: %s', requiredFields{iField});
    end
end

%% Validate matrix dimensions

assert(isequal(size(A.meanMatrix), [3 8]), ...
    'A.meanMatrix must be a 3x8 matrix.');

assert(isequal(size(A.maxMatrix), [3 8]), ...
    'A.maxMatrix must be a 3x8 matrix.');

assert(isequal(size(A.aucMatrix), [3 8]), ...
    'A.aucMatrix must be a 3x8 matrix.');

%% Convert condition-direction matrices to 24-element vectors
%
% The transpose followed by reshape preserves the ordering:
%
%   Cond1 Dir1-8
%   Cond2 Dir1-8
%   Cond3 Dir1-8
%

meanVals = reshape(A.meanMatrix.', 1, []);
maxVals  = reshape(A.maxMatrix.',  1, []);
aucVals  = reshape(A.aucMatrix.',  1, []);

%% Build output structure
Data = struct( ...
    'ICD_mean', num2cell(meanVals), ...
    'ICD_max',  num2cell(maxVals), ...
    'ICD_auc',  num2cell(aucVals));

end