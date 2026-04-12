function [f,J] = model_spm_fx_lfp_traj2d(x,u,P,M)
% model_spm_fx_lfp_traj2d  Hybrid DCM node model for two neural and two trajectory sources.
% FORMAT [f,J] = model_spm_fx_lfp_traj2d(x,u,P,M)
%
% Source order assumption
%   1. AreaS
%   2. AreaK
%   3. TrajS
%   4. TrajK
%
% Model idea
%   - AreaS and AreaK use the standard LFP neural-mass dynamics.
%   - TrajS and TrajK are latent trajectory nodes with simple second-order
%     dynamics embedded in the same 13-state source layout used by SPM.
%   - Crossed couplings follow the joint-action loop:
%         AreaS -> TrajS -> AreaK -> TrajK -> AreaS
%
% Notes
%   - This is a first hybrid implementation intended to make the modelling
%     assumptions explicit.
%   - The trajectory nodes only use states 1 and 4 as position and velocity;
%     all other states are softly damped to zero for numerical stability.
%   - Standard neural extrinsic coupling through A{1:3} is only applied
%     between the two neural sources.
%   - Extra trajectory couplings can be fixed in M.traj2d or, if desired in
%     the future, promoted to learnable parameters by defining priors on
%     P.Traj in a custom prior function.

try, P.H; catch, P.H = 0; end

x = spm_unvec(x,M.x);
nSources = size(x,1);
nStatesPerSource = size(x,2);

if nSources ~= 4
    error('%s expects exactly 4 sources ordered as [AreaS AreaK TrajS TrajK].', mfilename);
end

if nStatesPerSource < 13
    error('%s expects at least 13 states per source.', mfilename);
end

defaultExtrinsicRates = [32 16 4];
defaultIntrinsicRates = [1 1 1/2 1/2 1/32]*128;
defaultDelayConstants = [2 4];
defaultReceptorDensities = [8 32];
defaultSynapticConstants = [4 16];
defaultSigmoidParameters = [1 2];

if isfield(M,'pF')
    try, defaultExtrinsicRates = M.pF.E; end
    try, defaultIntrinsicRates = M.pF.H; end
    try, defaultDelayConstants = M.pF.D; end
    try, defaultReceptorDensities = M.pF.G; end
    try, defaultSynapticConstants = M.pF.T; end
    try, defaultSigmoidParameters = M.pF.R; end
end

A{1} = exp(P.A{1})*defaultExtrinsicRates(1);
A{2} = exp(P.A{2})*defaultExtrinsicRates(2);
A{3} = exp(P.A{3})*defaultExtrinsicRates(3);
C = exp(P.C);
intrinsicGainMatrix = exp(P.H)*diag(defaultIntrinsicRates);

timeConstantExcitatory = defaultSynapticConstants(1)/1000*exp(P.T(:,1));
timeConstantInhibitory = defaultSynapticConstants(2)/1000*exp(P.T(:,2));
timeConstantPotassium = 512/1000;
receptorDensityExcitatory = defaultReceptorDensities(1)*exp(P.G);
receptorDensityInhibitory = defaultReceptorDensities(2);

sigmoidParameters = defaultSigmoidParameters.*exp(P.R);
x = x';
shiftedState = x;
shiftedState(1,:) = shiftedState(1,:) - x(13,:);
firingRate = 1./(1 + exp(-sigmoidParameters(1)*(shiftedState - sigmoidParameters(2)))) ...
    - 1./(1 + exp(sigmoidParameters(1)*sigmoidParameters(2)));
firingRateDerivative = sigmoidParameters(1) ...
    * exp(-sigmoidParameters(1)*(max(shiftedState,-128) - sigmoidParameters(2))) ...
    ./ (1 + exp(-sigmoidParameters(1)*(shiftedState - sigmoidParameters(2)))).^2;

if isfield(M,'u')
    exogenousInput = u(:)*32;
else
    exogenousInput = C*u(:);
end

timeConstantInhibitory = 4/1000 + timeConstantInhibitory;

neuralSourceIndices = [1 2];
trajectorySourceIndices = [3 4];
areaToTrajectoryMap = [1 2];
trajectoryToAreaMap = [2 1];

trajectorySettings = getTrajectorySettings(M,P);

motionBySource = cell(1,nSources);
jacobianBlocks = cell(nSources,nSources);

for sourceIdx = neuralSourceIndices
    excitatoryTimeConstant = timeConstantExcitatory(sourceIdx);
    inhibitoryTimeConstant = timeConstantInhibitory(sourceIdx);
    excitatoryDecay = -2/excitatoryTimeConstant;
    inhibitoryDecay = -2/inhibitoryTimeConstant;
    excitatoryStiffness = -1/(excitatoryTimeConstant^2);
    inhibitoryStiffness = -1/(inhibitoryTimeConstant^2);

    localJacobian = [0   0   0   1   0   0   0   0   0   0   0   0   0
                     0   0   0   0   1   0   0   0   0   0   0   0   0
                     0   0   0   0   0   1   0   0   0   0   0   0   0
                     excitatoryStiffness  0   0   excitatoryDecay  0   0   0   0   0   0   0   0   0
                     0   excitatoryStiffness  0   0   excitatoryDecay  0   0   0   0   0   0   0   0
                     0   0   inhibitoryStiffness  0   0   inhibitoryDecay  0   0   0   0   0   0   0
                     0   0   0   0   0   0   0   1   0   0   0   0   0
                     0   0   0   0   0   0   excitatoryStiffness   excitatoryDecay  0   0   0   0   0
                     0   0   0   0   1   -1  0   0   0   0   0   0   0
                     0   0   0   0   0   0   0   0   0   0   1   0   0
                     0   0   0   0   0   0   0   0   0   inhibitoryStiffness   inhibitoryDecay  0   0
                     0   0   0   0   0   0   0   1   0   0  -1   0   0
                     0   0   0   0   0   0   0   0   0   0   0   0  -1/timeConstantPotassium];

    excitatoryGain = receptorDensityExcitatory(sourceIdx)/excitatoryTimeConstant;
    inhibitoryGain = receptorDensityInhibitory/inhibitoryTimeConstant;
    intrinsicDriveJacobian = [0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   excitatoryGain*intrinsicGainMatrix(sourceIdx,1)  0   0   0   0
                              excitatoryGain*intrinsicGainMatrix(sourceIdx,2)  0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   inhibitoryGain*intrinsicGainMatrix(sourceIdx,4)  0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   excitatoryGain*intrinsicGainMatrix(sourceIdx,3)  0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   inhibitoryGain*intrinsicGainMatrix(sourceIdx,5)   0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              4/timeConstantPotassium 0 0 0 0 0 0 0 0 0 0 0 0];

    stateDerivativeJacobian = diag(firingRateDerivative(:,sourceIdx));
    stateDerivativeJacobian(1,13) = -stateDerivativeJacobian(1,1);
    inputJacobian = sparse(4,1,excitatoryGain,nStatesPerSource,1);

    motionBySource{sourceIdx} = localJacobian*x(:,sourceIdx) ...
        + intrinsicDriveJacobian*firingRate(:,sourceIdx) ...
        + inputJacobian*exogenousInput(sourceIdx);
    jacobianBlocks{sourceIdx,sourceIdx} = localJacobian + intrinsicDriveJacobian*stateDerivativeJacobian;

    for afferentIdx = neuralSourceIndices
        if afferentIdx == sourceIdx
            continue;
        end

        forwardInputGain = excitatoryGain*(A{1}(sourceIdx,afferentIdx) + A{3}(sourceIdx,afferentIdx));
        backwardInputGain = excitatoryGain*(A{2}(sourceIdx,afferentIdx) + A{3}(sourceIdx,afferentIdx));
        inhibitoryInputGain = inhibitoryGain*(A{2}(sourceIdx,afferentIdx) + A{3}(sourceIdx,afferentIdx));

        afferentDriveJacobian = [0   0   0   0   0   0   0   0   0   0   0   0   0
                                 0   0   0   0   0   0   0   0   0   0   0   0   0
                                 0   0   0   0   0   0   0   0   0   0   0   0   0
                                 0   0   0   0   0   0   0   0   forwardInputGain   0   0   0   0
                                 0   0   0   0   0   0   0   0   backwardInputGain  0   0   0   0
                                 0   0   0   0   0   0   0   0   0   0   0   0   0
                                 0   0   0   0   0   0   0   0   0   0   0   0   0
                                 0   0   0   0   0   0   0   0   inhibitoryInputGain 0  0   0   0
                                 0   0   0   0   0   0   0   0   0   0   0   0   0
                                 0   0   0   0   0   0   0   0   0   0   0   0   0
                                 0   0   0   0   0   0   0   0   0   0   0   0   0
                                 0   0   0   0   0   0   0   0   0   0   0   0   0
                                 0   0   0   0   0   0   0   0   0   0   0   0   0];

        motionBySource{sourceIdx} = motionBySource{sourceIdx} + afferentDriveJacobian*firingRate(:,afferentIdx);
        jacobianBlocks{sourceIdx,afferentIdx} = afferentDriveJacobian*diag(firingRateDerivative(:,afferentIdx));
    end
end

for trajectoryNodeOffset = 1:numel(trajectorySourceIndices)
    trajectorySourceIdx = trajectorySourceIndices(trajectoryNodeOffset);
    drivingAreaIdx = areaToTrajectoryMap(trajectoryNodeOffset);

    trajectoryFrequency = trajectorySettings.angularFrequency(trajectoryNodeOffset);
    trajectoryDamping = trajectorySettings.dampingRatio(trajectoryNodeOffset);
    areaToTrajectoryGain = trajectorySettings.areaToTrajectoryGain(trajectoryNodeOffset);
    unusedStateDecay = trajectorySettings.unusedStateDecay;

    trajectoryMotion = zeros(nStatesPerSource,1);
    trajectoryJacobianSelf = -unusedStateDecay*speye(nStatesPerSource);
    trajectoryJacobianSelf(1,1) = 0;
    trajectoryJacobianSelf(1,4) = 1;
    trajectoryJacobianSelf(4,1) = -(trajectoryFrequency^2);
    trajectoryJacobianSelf(4,4) = -2*trajectoryDamping*trajectoryFrequency;

    trajectoryMotion(1) = x(4,trajectorySourceIdx);
    trajectoryMotion(4) = -(trajectoryFrequency^2)*x(1,trajectorySourceIdx) ...
        - 2*trajectoryDamping*trajectoryFrequency*x(4,trajectorySourceIdx) ...
        + areaToTrajectoryGain*firingRate(9,drivingAreaIdx);

    trajectoryMotion(2:end) = trajectoryMotion(2:end) - unusedStateDecay*x(2:end,trajectorySourceIdx);
    trajectoryMotion(4) = trajectoryMotion(4) + unusedStateDecay*x(4,trajectorySourceIdx);

    motionBySource{trajectorySourceIdx} = trajectoryMotion;
    jacobianBlocks{trajectorySourceIdx,trajectorySourceIdx} = trajectoryJacobianSelf;

    sourceToTrajectoryJacobian = zeros(nStatesPerSource,nStatesPerSource);
    sourceToTrajectoryJacobian(4,9) = areaToTrajectoryGain*firingRateDerivative(9,drivingAreaIdx);
    jacobianBlocks{trajectorySourceIdx,drivingAreaIdx} = sourceToTrajectoryJacobian;
end

for trajectoryNodeOffset = 1:numel(trajectorySourceIndices)
    sourceTrajectoryIdx = trajectorySourceIndices(trajectoryNodeOffset);
    targetAreaIdx = trajectoryToAreaMap(trajectoryNodeOffset);

    targetExcitatoryGain = receptorDensityExcitatory(targetAreaIdx)/timeConstantExcitatory(targetAreaIdx);
    trajectoryToAreaGain = trajectorySettings.trajectoryToAreaGain(trajectoryNodeOffset);

    motionBySource{targetAreaIdx}(4) = motionBySource{targetAreaIdx}(4) ...
        + targetExcitatoryGain*trajectoryToAreaGain*x(1,sourceTrajectoryIdx);

    if isempty(jacobianBlocks{targetAreaIdx,sourceTrajectoryIdx})
        jacobianBlocks{targetAreaIdx,sourceTrajectoryIdx} = zeros(nStatesPerSource,nStatesPerSource);
    end
    jacobianBlocks{targetAreaIdx,sourceTrajectoryIdx}(4,1) = ...
        jacobianBlocks{targetAreaIdx,sourceTrajectoryIdx}(4,1) ...
        + targetExcitatoryGain*trajectoryToAreaGain;
end

f = zeros(nSources*nStatesPerSource,1);
fullJacobian = sparse(nSources*nStatesPerSource,nSources*nStatesPerSource);

for sourceIdx = 1:nSources
    rowIndices = (1:nSources:nStatesPerSource*nSources) + (sourceIdx - 1);
    f(rowIndices,1) = motionBySource{sourceIdx};
    for afferentIdx = 1:nSources
        columnIndices = (1:nSources:nStatesPerSource*nSources) + (afferentIdx - 1);
        if ~isempty(jacobianBlocks{sourceIdx,afferentIdx})
            fullJacobian(rowIndices,columnIndices) = jacobianBlocks{sourceIdx,afferentIdx};
        end
    end
end

extrinsicDelayMatrix = defaultDelayConstants(2).*exp(P.D)/1000;
intrinsicDelayMatrix = defaultDelayConstants(1).*exp(P.I)/1000;
extrinsicDelayMatrix = (eye(nSources,nSources) - 1).*extrinsicDelayMatrix;
intrinsicDelayMatrix = (eye(nStatesPerSource,nStatesPerSource) - 1)*intrinsicDelayMatrix;
extrinsicDelayMatrix = kron(ones(nStatesPerSource,nStatesPerSource),extrinsicDelayMatrix);
intrinsicDelayMatrix = kron(intrinsicDelayMatrix,eye(nSources,nSources));
delayOperator = intrinsicDelayMatrix + extrinsicDelayMatrix;

delayOperator = spm_inv(speye(nSources*nStatesPerSource,nSources*nStatesPerSource) - delayOperator.*fullJacobian);
f = delayOperator*f;
J = delayOperator*fullJacobian;
end


function trajectorySettings = getTrajectorySettings(M,P)
trajectorySettings = struct();
trajectorySettings.angularFrequency = [4 4];
trajectorySettings.dampingRatio = [1 1];
trajectorySettings.areaToTrajectoryGain = [1 1];
trajectorySettings.trajectoryToAreaGain = [1 1];
trajectorySettings.unusedStateDecay = 32;

if isfield(M,'traj2d')
    if isfield(M.traj2d,'angularFrequency')
        trajectorySettings.angularFrequency = M.traj2d.angularFrequency;
    end
    if isfield(M.traj2d,'dampingRatio')
        trajectorySettings.dampingRatio = M.traj2d.dampingRatio;
    end
    if isfield(M.traj2d,'areaToTrajectoryGain')
        trajectorySettings.areaToTrajectoryGain = M.traj2d.areaToTrajectoryGain;
    end
    if isfield(M.traj2d,'trajectoryToAreaGain')
        trajectorySettings.trajectoryToAreaGain = M.traj2d.trajectoryToAreaGain;
    end
    if isfield(M.traj2d,'unusedStateDecay')
        trajectorySettings.unusedStateDecay = M.traj2d.unusedStateDecay;
    end
end

if isfield(P,'Traj')
    if isfield(P.Traj,'angularFrequency')
        trajectorySettings.angularFrequency = exp(P.Traj.angularFrequency(:))';
    end
    if isfield(P.Traj,'dampingRatio')
        trajectorySettings.dampingRatio = exp(P.Traj.dampingRatio(:))';
    end
    if isfield(P.Traj,'areaToTrajectoryGain')
        trajectorySettings.areaToTrajectoryGain = exp(P.Traj.areaToTrajectoryGain(:))';
    end
    if isfield(P.Traj,'trajectoryToAreaGain')
        trajectorySettings.trajectoryToAreaGain = exp(P.Traj.trajectoryToAreaGain(:))';
    end
end
end
