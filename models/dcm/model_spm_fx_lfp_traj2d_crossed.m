function [stateDerivative, jacobianMatrix] = model_spm_fx_lfp_traj2d_crossed(stateVector, inputVector, parameterStruct, modelStruct)
% model_spm_fx_lfp_traj2d_crossed
% Hybrid 4-source model with two neural areas and two latent trajectory nodes.
%
% Source order
%   1. AreaS
%   2. AreaK
%   3. TrajS
%   4. TrajK
%
% Architecture
%   Only A{1} is used. A{2} and A{3} are ignored by design.
%
%   Structural loop encoded in A{1}:
%     A{1}(1,4) : TrajK -> AreaS
%     A{1}(2,3) : TrajS -> AreaK
%     A{1}(3,1) : AreaS -> TrajS
%     A{1}(4,2) : AreaK -> TrajK
%
%   Modulations are expected to be folded into the effective P.A{1} by
%   spm_gen_Q before this state equation is called. This makes it natural
%   to let B blocks modulate only the crossed trajectory-to-area influences.

try, parameterStruct.H; catch, parameterStruct.H = 0; end

stateBySource = spm_unvec(stateVector, modelStruct.x);
nSources = size(stateBySource, 1);
nStatesPerSource = size(stateBySource, 2);

if nSources ~= 4
    error('%s expects 4 sources ordered as [AreaS AreaK TrajS TrajK].', mfilename);
end

if nStatesPerSource < 13
    error('%s expects at least 13 states per source.', mfilename);
end

defaultExtrinsicRates = [32 16 4];
defaultIntrinsicRates = [1 1 1/2 1/2 1/32] * 128;
defaultDelayConstants = [2 4];
defaultReceptorDensities = [8 32];
defaultSynapticConstants = [4 16];
defaultSigmoidParameters = [1 2];

if isfield(modelStruct, 'pF')
    try, defaultExtrinsicRates = modelStruct.pF.E; end
    try, defaultIntrinsicRates = modelStruct.pF.H; end
    try, defaultDelayConstants = modelStruct.pF.D; end
    try, defaultReceptorDensities = modelStruct.pF.G; end
    try, defaultSynapticConstants = modelStruct.pF.T; end
    try, defaultSigmoidParameters = modelStruct.pF.R; end
end

forwardMask = local_get_forward_mask(modelStruct, nSources);
effectiveForwardCoupling = forwardMask .* exp(parameterStruct.A{1});

intrinsicGainMatrix = exp(parameterStruct.H) * diag(defaultIntrinsicRates);
timeConstantExcitatory = defaultSynapticConstants(1) / 1000 * exp(parameterStruct.T(:,1));
timeConstantInhibitory = defaultSynapticConstants(2) / 1000 * exp(parameterStruct.T(:,2));
timeConstantPotassium = 512 / 1000;
receptorDensityExcitatory = defaultReceptorDensities(1) * exp(parameterStruct.G);
receptorDensityInhibitory = defaultReceptorDensities(2);

sigmoidParameters = defaultSigmoidParameters .* exp(parameterStruct.R);
stateByPopulation = stateBySource';
shiftedState = stateByPopulation;
shiftedState(1,:) = shiftedState(1,:) - stateByPopulation(13,:);
firingRate = 1 ./ (1 + exp(-sigmoidParameters(1) * (shiftedState - sigmoidParameters(2)))) ...
    - 1 ./ (1 + exp(sigmoidParameters(1) * sigmoidParameters(2)));
firingRateDerivative = sigmoidParameters(1) ...
    * exp(-sigmoidParameters(1) * (max(shiftedState, -128) - sigmoidParameters(2))) ...
    ./ (1 + exp(-sigmoidParameters(1) * (shiftedState - sigmoidParameters(2)))).^2;

exogenousInput = local_resolve_exogenous_input(inputVector, parameterStruct, modelStruct, nSources);
timeConstantInhibitory = 4 / 1000 + timeConstantInhibitory;
trajectorySettings = local_get_trajectory_settings(modelStruct);

neuralSourceIndices = [1 2];
motionBySource = cell(1, nSources);
jacobianBlocks = cell(nSources, nSources);

for neuralSourceIndex = neuralSourceIndices
    excitatoryTimeConstant = timeConstantExcitatory(neuralSourceIndex);
    inhibitoryTimeConstant = timeConstantInhibitory(neuralSourceIndex);
    excitatoryDecay = -2 / excitatoryTimeConstant;
    inhibitoryDecay = -2 / inhibitoryTimeConstant;
    excitatoryStiffness = -1 / (excitatoryTimeConstant ^ 2);
    inhibitoryStiffness = -1 / (inhibitoryTimeConstant ^ 2);

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
                     0   0   0   0   0   0   0   0   0   0   0   0  -1 / timeConstantPotassium];

    excitatoryGain = receptorDensityExcitatory(neuralSourceIndex) / excitatoryTimeConstant;
    inhibitoryGain = receptorDensityInhibitory / inhibitoryTimeConstant;
    intrinsicDriveJacobian = [0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   excitatoryGain * intrinsicGainMatrix(neuralSourceIndex,1)  0   0   0   0
                              excitatoryGain * intrinsicGainMatrix(neuralSourceIndex,2)  0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   inhibitoryGain * intrinsicGainMatrix(neuralSourceIndex,4)  0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   excitatoryGain * intrinsicGainMatrix(neuralSourceIndex,3)  0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              0   0   0   0   0   0   0   0   0   0   0   inhibitoryGain * intrinsicGainMatrix(neuralSourceIndex,5)   0
                              0   0   0   0   0   0   0   0   0   0   0   0   0
                              4 / timeConstantPotassium 0 0 0 0 0 0 0 0 0 0 0 0];

    stateDerivativeJacobian = diag(firingRateDerivative(:,neuralSourceIndex));
    stateDerivativeJacobian(1,13) = -stateDerivativeJacobian(1,1);
    inputJacobian = sparse(4,1,excitatoryGain,nStatesPerSource,1);

    motionBySource{neuralSourceIndex} = localJacobian * stateByPopulation(:,neuralSourceIndex) ...
        + intrinsicDriveJacobian * firingRate(:,neuralSourceIndex) ...
        + inputJacobian * exogenousInput(neuralSourceIndex);
    jacobianBlocks{neuralSourceIndex, neuralSourceIndex} = localJacobian + intrinsicDriveJacobian * stateDerivativeJacobian;
end

crossedTrajectoryToAreaLinks = [1 4; 2 3];
for linkIndex = 1:size(crossedTrajectoryToAreaLinks, 1)
    targetAreaIndex = crossedTrajectoryToAreaLinks(linkIndex, 1);
    sourceTrajectoryIndex = crossedTrajectoryToAreaLinks(linkIndex, 2);
    targetExcitatoryGain = receptorDensityExcitatory(targetAreaIndex) / timeConstantExcitatory(targetAreaIndex);
    trajectoryCouplingGain = trajectorySettings.trajectoryToAreaBaseGain(linkIndex) ...
        * effectiveForwardCoupling(targetAreaIndex, sourceTrajectoryIndex);

    motionBySource{targetAreaIndex}(4) = motionBySource{targetAreaIndex}(4) ...
        + targetExcitatoryGain * trajectoryCouplingGain * stateByPopulation(1, sourceTrajectoryIndex);

    trajectoryToAreaJacobian = zeros(nStatesPerSource, nStatesPerSource);
    trajectoryToAreaJacobian(4,1) = targetExcitatoryGain * trajectoryCouplingGain;
    jacobianBlocks{targetAreaIndex, sourceTrajectoryIndex} = trajectoryToAreaJacobian;
end

areaToTrajectoryLinks = [3 1; 4 2];
for linkIndex = 1:size(areaToTrajectoryLinks, 1)
    targetTrajectoryIndex = areaToTrajectoryLinks(linkIndex, 1);
    sourceAreaIndex = areaToTrajectoryLinks(linkIndex, 2);

    trajectoryFrequency = trajectorySettings.angularFrequency(linkIndex);
    trajectoryDamping = trajectorySettings.dampingRatio(linkIndex);
    areaToTrajectoryGain = trajectorySettings.areaToTrajectoryBaseGain(linkIndex) ...
        * effectiveForwardCoupling(targetTrajectoryIndex, sourceAreaIndex);
    unusedStateDecay = trajectorySettings.unusedStateDecay;

    trajectoryMotion = zeros(nStatesPerSource, 1);
    trajectoryJacobianSelf = -unusedStateDecay * speye(nStatesPerSource);
    trajectoryJacobianSelf(1,1) = 0;
    trajectoryJacobianSelf(1,4) = 1;
    trajectoryJacobianSelf(4,1) = -(trajectoryFrequency ^ 2);
    trajectoryJacobianSelf(4,4) = -2 * trajectoryDamping * trajectoryFrequency;

    trajectoryMotion(1) = stateByPopulation(4, targetTrajectoryIndex);
    trajectoryMotion(4) = -(trajectoryFrequency ^ 2) * stateByPopulation(1, targetTrajectoryIndex) ...
        - 2 * trajectoryDamping * trajectoryFrequency * stateByPopulation(4, targetTrajectoryIndex) ...
        + areaToTrajectoryGain * firingRate(9, sourceAreaIndex);

    trajectoryMotion(2:end) = trajectoryMotion(2:end) - unusedStateDecay * stateByPopulation(2:end, targetTrajectoryIndex);
    trajectoryMotion(4) = trajectoryMotion(4) + unusedStateDecay * stateByPopulation(4, targetTrajectoryIndex);

    motionBySource{targetTrajectoryIndex} = trajectoryMotion;
    jacobianBlocks{targetTrajectoryIndex, targetTrajectoryIndex} = trajectoryJacobianSelf;

    areaToTrajectoryJacobian = zeros(nStatesPerSource, nStatesPerSource);
    areaToTrajectoryJacobian(4,9) = areaToTrajectoryGain * firingRateDerivative(9, sourceAreaIndex);
    jacobianBlocks{targetTrajectoryIndex, sourceAreaIndex} = areaToTrajectoryJacobian;
end

stateDerivative = zeros(nSources * nStatesPerSource, 1);
fullJacobian = sparse(nSources * nStatesPerSource, nSources * nStatesPerSource);

for sourceIndex = 1:nSources
    rowIndices = (1:nSources:nStatesPerSource * nSources) + (sourceIndex - 1);
    stateDerivative(rowIndices, 1) = motionBySource{sourceIndex};
    for afferentSourceIndex = 1:nSources
        columnIndices = (1:nSources:nStatesPerSource * nSources) + (afferentSourceIndex - 1);
        if ~isempty(jacobianBlocks{sourceIndex, afferentSourceIndex})
            fullJacobian(rowIndices, columnIndices) = jacobianBlocks{sourceIndex, afferentSourceIndex};
        end
    end
end

extrinsicDelayMatrix = defaultDelayConstants(2) .* exp(parameterStruct.D) / 1000;
intrinsicDelayMatrix = defaultDelayConstants(1) .* exp(parameterStruct.I) / 1000;
extrinsicDelayMatrix = (eye(nSources, nSources) - 1) .* extrinsicDelayMatrix;
intrinsicDelayMatrix = (eye(nStatesPerSource, nStatesPerSource) - 1) * intrinsicDelayMatrix;
extrinsicDelayMatrix = kron(ones(nStatesPerSource, nStatesPerSource), extrinsicDelayMatrix);
intrinsicDelayMatrix = kron(intrinsicDelayMatrix, eye(nSources, nSources));
delayOperator = intrinsicDelayMatrix + extrinsicDelayMatrix;

delayOperator = spm_inv(speye(nSources * nStatesPerSource, nSources * nStatesPerSource) - delayOperator .* fullJacobian);
stateDerivative = delayOperator * stateDerivative;
jacobianMatrix = delayOperator * fullJacobian;
end


function forwardMask = local_get_forward_mask(modelStruct, nSources)
forwardMask = zeros(nSources, nSources);

if isfield(modelStruct, 'architecture') && isfield(modelStruct.architecture, 'A') ...
        && numel(modelStruct.architecture.A) >= 1 && ~isempty(modelStruct.architecture.A{1})
    forwardMask = double(modelStruct.architecture.A{1} ~= 0);
end
end


function exogenousInput = local_resolve_exogenous_input(inputVector, parameterStruct, modelStruct, nSources)
if nargin < 1 || isempty(inputVector)
    exogenousInput = zeros(nSources, 1);
    return;
end

if isfield(modelStruct, 'u') && ~isempty(modelStruct.u) && spm_length(modelStruct.u) > 0
    exogenousInput = inputVector(:) * 32;
    return;
end

if isfield(parameterStruct, 'C') && ~isempty(parameterStruct.C)
    exogenousInput = exp(parameterStruct.C) * inputVector(:);
    return;
end

exogenousInput = zeros(nSources, 1);
end


function trajectorySettings = local_get_trajectory_settings(modelStruct)
trajectorySettings = struct();
trajectorySettings.angularFrequency = [5 5];
trajectorySettings.dampingRatio = [1.0 1.0];
trajectorySettings.areaToTrajectoryBaseGain = [1.6 1.6];
trajectorySettings.trajectoryToAreaBaseGain = [0.5 0.5];
trajectorySettings.unusedStateDecay = 32;

if isfield(modelStruct, 'crossedTraj2d')
    if isfield(modelStruct.crossedTraj2d, 'angularFrequency')
        trajectorySettings.angularFrequency = modelStruct.crossedTraj2d.angularFrequency;
    end
    if isfield(modelStruct.crossedTraj2d, 'dampingRatio')
        trajectorySettings.dampingRatio = modelStruct.crossedTraj2d.dampingRatio;
    end
    if isfield(modelStruct.crossedTraj2d, 'areaToTrajectoryBaseGain')
        trajectorySettings.areaToTrajectoryBaseGain = modelStruct.crossedTraj2d.areaToTrajectoryBaseGain;
    end
    if isfield(modelStruct.crossedTraj2d, 'trajectoryToAreaBaseGain')
        trajectorySettings.trajectoryToAreaBaseGain = modelStruct.crossedTraj2d.trajectoryToAreaBaseGain;
    end
    if isfield(modelStruct.crossedTraj2d, 'unusedStateDecay')
        trajectorySettings.unusedStateDecay = modelStruct.crossedTraj2d.unusedStateDecay;
    end
end
end
