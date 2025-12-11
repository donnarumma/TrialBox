function gx = spm_gx_multiERP(x, u, P, M)
% spm_gx_multierp  Observation function for multi-channel, multi-source LFP DCM
%
% FORMAT gx = spm_gx_multiERP(x, u, P, M)
%
% x : [nsources x nstates] neural states
% u : inputs (unused here)
% P : parameter structure with fields:
%     .J - readout weights for neural states (1 x nstates)
% M : model structure with field:
%     .L - channel-to-source mapping matrix (nchannels x nsources)
%     .x - baseline state [nsources x nstates]
%
% OUTPUT:
% gx : predicted LFP/EEG/MEG signal at each channel

% Ensure J is a column vector
J                   = P.J(:);                   % [nstates x 1]

% Project each source's state onto observation dimension
% Result: [nsources x 1]
dx                  = x - M.x;                  % deviation from baseline
source_signals      = dx * J;                   % project each source to 1D

% Map source signals to channels using lead field matrix M.L
% M.L is [nchannels x nsources], source_signals is [nsources x 1]
gx                  = P.L * source_signals;     % [nchannels x 1]
end
