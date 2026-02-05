% =========================================================================
% DATA STRUCTURE REQUIREMENTS (MANDATORY)
%
% This pipeline ASSUMES that each element of the input subject structures
% (e.g. X_struct, Y_struct) MUST contain
% the following fields with EXACT NAMES and SEMANTICS:
% Each structure is a collection of trials of aggregated sessions
%
% REQUIRED FIELDS (case-sensitive):
%
%   - trialTypeDir   : scalar numeric
%       Direction label of the trial (e.g. 1..max_dir)
%
%   - trialTypeCond  : scalar numeric
%       Condition label of the trial (e.g. 1..max_cond)
%
%   - Spikes         : numeric matrix [N_neurons x N_timebins per trial]
%       Neural activity for the trial.
%       ALL trials used in this analysis MUST have the SAME time length.
%
%   - Manifold       : numeric matrix [D x N_timebins, per trial]
%       Low-dimensional representation (e.g. CEBRA embedding),
%       temporally aligned with Spikes.
%
% IMPORTANT NOTES:
% - Field names MUST match exactly (no aliases, no legacy names).
% - Fields must NOT be empty.
% - trialTypeDir and trialTypeCond MUST be scalar values.
% - Trials violating these assumptions will cause UNDEFINED BEHAVIOR
%   or runtime errors in downstream functions (e.g. filter_data).
%
% If your raw data uses different field names (e.g. trial_type_d,
% trialTypeCon, etc.), YOU MUST NORMALIZE THEM *BEFORE* calling this code.
%
% This file DOES NOT perform automatic field-name correction by design.
% =========================================================================
