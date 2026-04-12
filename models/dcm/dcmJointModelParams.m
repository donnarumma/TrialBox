function par = dcmJointModelParams(parSource)
% dcmJointModelParams  Default parameters for dcmJointModel.
%
% PURPOSE
%   Provide reusable defaults for session-level DCM fitting of the joint
%   monkey model before the call to csdModel and the variational inversion.
%
% PARAMETERS
%   mstep                 : M-step size used by the inversion routine.
%   isdebug               : debug flag passed to csdModel.
%   Nmax                  : maximum number of frequency bins used internally.
%   whichmodel            : fallback model index when A/B/C are not given.
%   custom_model          : optional custom prior function.
%   InField               : trial field used as spectral input, usually CSD.
%   donlfp                : whether to use the custom LFP priors.
%   isdemo                : selected direction-condition indices.
%   S_PRIORS              : reserved hook for monkey S priors.
%   K_PRIORS              : reserved hook for monkey K priors.
%   Hz                    : frequency axis used by the model.
%   A                     : baseline architecture matrices.
%   B                     : modulatory architecture matrices.
%   C                     : driving-input architecture matrix.
%   L                     : observation gain matrix for source/channel mismatch.
%   description           : model label saved in the fitted DCM.
%   stateEquationFunction : state equation assigned to M.f.
%   traj2d                : optional trajectory settings for trajectory models.
%   crossedTraj2d         : optional crossed trajectory settings for ModelM4.

par.mstep = 8;
par.isdebug = false;
par.Nmax = 128;
par.whichmodel = 5;
par.custom_model = [];
par.InField = 'CSD';
par.donlfp = false;
par.isdemo = 1:24;
par.S_PRIORS = [];
par.K_PRIORS = [];
par.Hz = (1:64)';
par.A = [];
par.B = [];
par.C = [];
par.L = [];
par.description = '';
par.stateEquationFunction = 'model_spm_fx_lfp';
par.traj2d = [];
par.crossedTraj2d = [];

if nargin < 1 || isempty(parSource)
    return;
end

fieldNames = fieldnames(parSource);
for fieldIndex = 1:numel(fieldNames)
    par.(fieldNames{fieldIndex}) = parSource.(fieldNames{fieldIndex});
end
end
