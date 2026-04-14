function   M_CSD=csdModel(par)%nConditions,nSources,nChannels,whichmodel,customMode,donlfp)
nConditions = par.nConditions;
nSources    = par.nSources;
nChannels   = par.nChannels; 
whichmodel  = par.whichmodel; 
customMode  = par.customMode;
donlfp      = par.donlfp;  
mstep       = par.mstep;
isdebug     = par.isdebug;
Nmax        = par.Nmax;
stateEquationFunction = getOption(par,'stateEquationFunction','model_spm_fx_lfp');
% function M_CSD=csdModel(nConditions,nSources,nChannels,whichmodel,customMode,donlfp)
% Ntrials: number of trials

% number of sources and LFP channels (usually the same)
%--------------------------------------------------------------------------
try
    nSo   = nSources;
catch
    nSo   = 2;             % number of sources
end
try
    nCh  = nChannels; 
catch
    nCh  = nSo; %2; % number of channels
end
% specify network (connections)
%--------------------------------------------------------------------------
try 
    A=par.A;
    B=par.B;
    C=par.C;
catch
    switch whichmodel
        case 1
            A{1}  = tril(ones(nSo,nSo),-1);       % a forward connection
            A{2}  = triu(ones(nSo,nSo),+1);       % a backward connection
            A{3}  = sparse(nSo,nSo);              % lateral connections
            M.description='OM1';
        case 2
            % iss inihibitory spiny stellate 
            % ip  inhibitory piramidals
            % ii  inhibitory interneurons 
            % 
            A{1}  = ones(nSo,nSo)-eye(nSo);       % iss <-> iss
            A{2}  = ones(nSo,nSo)-eye(nSo);       %  ip <->  ip + ii <-> ii
            A{3}  = sparse(nSo,nSo);              % iss <-> iss + ip <-> ip + ii <-> ii -- empty 
            M.description='OM2';
        case 3
            % only A{3} 
            A{1}  = sparse(nSo,nSo);
            A{2}  = sparse(nSo,nSo);
            A{3}  = ones(nSo,nSo)-eye(nSo);
            M.description='OM3';
        case 4
            % only A{3} with autoconnections                
            A{1}  = sparse(nSo,nSo);
            A{2}  = sparse(nSo,nSo);
            A{3}  = ones(nSo,nSo);
            M.description='OM4';
        case 5 
            % only B interactions
            A{1}  = sparse(nSo,nSo);
            A{2}  = sparse(nSo,nSo);
            A{3}  = sparse(nSo,nSo);
            M.description='OM5';
        case 6
            % autoconnections A{3}
            A{1}  = sparse(nSo,nSo);
            A{2}  = sparse(nSo,nSo);
            A{3}  = eye(nSo,nSo);
            M.description='OM6';
        case 7
            % autoconnections A + B
            A{1}  = eye(nSo,nSo);
            A{2}  = eye(nSo,nSo);
            A{3}  = eye(nSo,nSo);
            M.description='OM7'; 
        case 8
            % connections A
            A{1}  = ones(nSo,nSo);
            A{2}  = ones(nSo,nSo);
            A{3}  = ones(nSo,nSo);
            M.description='OM8';            
        case 9
            % all autoconnections
            A{1}  = ones(nSo,nSo)-eye(nSo);
            A{2}  = ones(nSo,nSo)-eye(nSo);
            A{3}  = ones(nSo,nSo)-eye(nSo);
            M.description='OM9';
        case 10
            % autoconnections A + interactions B
            A{1}  = eye(nSo,nSo);
            A{2}  = eye(nSo,nSo);
            A{3}  = eye(nSo,nSo);
            M.description='OM10';
    end
    
    try
        nCo    = nConditions;
    catch
        nCo    = 2;                             % number of conditions
    end
    B           = cell(1,nCo);                  % trial-specific modulation -> see spm_gen_Q
    for iCo=1:nCo
        if whichmodel==7
            B{iCo}  = eye(nSo);
        else
            B{iCo}  = ones(nSo,nSo)-eye(nSo); 
        end
    end
    if nCo==1  %% WARNING if more trials of a unique type, B is unused
        B{1}=sparse(nSo,nSo);
    end
    C           = speye(nSo,nSo);               % sources receiving innovations
end
% get priors
%--------------------------------------------------------------------------
if ~donlfp 
    [pE,pC]     = spm_lfp_priors(A,B,C);        % neuronal priors
else
    [pE,pC]     = donnarumma_spm_lfp_priors(A,B,C);        % neuronal priors
end
[pE,pC]     = spm_ssr_priors(pE,pC);        % spectral priors
dipfitSpec  = struct('type','LFP','model','LFP','Ns',nSo,'Nc',nCh);
[pE,pC]     = spm_L_priors(dipfitSpec,pE,pC);      % spatial  priors
if ~isempty(customMode)
    [pE,pC]     = customMode(pE,pC,par);    % load custom priors
end
% Suppress channel noise
%--------------------------------------------------------------------------
pE.b        = pE.b - 16;
pE.c        = pE.c - 16;
% create LFP model
%--------------------------------------------------------------------------
% Standard DCM/CSD pipeline in SPM:
%   1) DONNARUMMA_spm_nlsi_GN composes M.IS and M.FS as
%      predicted_features = M.FS(M.IS(P,M,U), M)
%   2) M.IS = spm_csd_mtf builds condition-specific CSD predictions:
%      - spm_gen_Q applies between-condition effects in U.X to P.B
%        and therefore expects one P.B{i} per column of U.X
%      - spm_dcm_mtf linearises the neural model around steady state
%      - channel noise is then added in sensor space
%   3) M.f = model_spm_fx_lfp defines the neuronal state equations
%   4) M.g = spm_gx_erp maps hidden states to sensor space
%   5) M.FS = spm_fs_csd converts the predicted CSD into the feature vector
%      actually used by the variational inversion
%
% These hooks are standard in the DCM-for-CSD framework; we collect them in
% one place so the model definition is easier to read and customise.
% To use the hybrid neural/trajectory model, pass for example:
%   par.stateEquationFunction = 'model_spm_fx_lfp_traj2d';
%   par.traj2d = struct(...);
% or for the crossed four-source variant:
%   par.stateEquationFunction = 'model_spm_fx_lfp_traj2d_crossed';
%   par.crossedTraj2d = struct(...);
dcmCsdHooks = getStandardDcmCsdHooks();
dcmCsdHooks.stateEquationFunction = stateEquationFunction;

M.dipfit = dipfitSpec;
M.architecture = struct();
M.architecture.A = A;
M.architecture.B = B;
M.architecture.C = C;

M.IS        = dcmCsdHooks.predictionFunction;        % build predicted CSDs
M.FS        = dcmCsdHooks.featureSelectionFunction;  % convert CSDs into inversion features
M.g         = dcmCsdHooks.observerFunction;          % map hidden states to sensor space
%% alternative neural models
%% see -> spm_fx_cmc_tfm_gen    ( 8 neurons)  A{4} see slides - spm12/toolbox/NVC
%% see -> spm_fx_lfp            (13 neurons)
%% see -> spm_fx_fmri           ( 8 cells - BOLD)
%% see -> spm_fx_erp            ( 9 neurons)
M.f         = dcmCsdHooks.stateEquationFunction;     % define neuronal state equations
M.x         = sparse(nSo,13);     
M.n         = nSo*13;
M.pE        = pE;
M.pC        = pC;
M.m         = nSo;            % areas
M.l         = nCh;            % electrodes        
M.Hz        = par.Hz;         % (1:64)';

M.Nmax      = 256;
M.nograph   = 1;              % no plot

% case nChannels > nSources
if nSources<nChannels
    M.pE.L  =par.L;
    M.pC.L  =par.L*mean(pC.L);
    M.g     = 'spm_gx_multiERP';
end

M.mstep     = mstep;
M.isdebug   = isdebug;
M.Nmax      = Nmax;
M.dcmCsdHooks = dcmCsdHooks;
if isfield(par,'traj2d')
    M.traj2d = par.traj2d;
end
if isfield(par,'crossedTraj2d')
    M.crossedTraj2d = par.crossedTraj2d;
end

M_CSD       = M;


function dcmCsdHooks = getStandardDcmCsdHooks()
dcmCsdHooks = struct();
dcmCsdHooks.predictionFunction          = 'spm_csd_mtf';        % M.IS - build predicted CSDs
dcmCsdHooks.featureSelectionFunction    = 'spm_fs_csd';         % M.FS - convert predicted CSDs into inversion features
dcmCsdHooks.observerFunction            = 'spm_gx_erp';         % M.g  - map hidden states to sensor space
dcmCsdHooks.stateEquationFunction       = 'model_spm_fx_lfp';   % M.f  - define neuronal state equations
end


function optionValue = getOption(optionStruct, optionName, defaultValue)
if isfield(optionStruct, optionName)
    optionValue = optionStruct.(optionName);
else
    optionValue = defaultValue;
end
end

end
