function   M_CSD=csdModel(par)%nConditions,nSources,nChannels,whichmodel,customMode,donlfp)
nConditions = par.nConditions;
nSources    = par.nSources;
nChannels   = par.nChannels; 
whichmodel  = par.whichmodel; 
customMode  = par.customMode;
donlfp      = par.donlfp;  

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
        case 2
            % iss inihibitory spiny stellate 
            % ip  inhibitory piramidals
            % ii  inhibitory interneurons 
            % 
            A{1}  = ones(nSo,nSo)-eye(nSo);       % iss <-> iss
            A{2}  = ones(nSo,nSo)-eye(nSo);       %  ip <->  ip + ii <-> ii
            A{3}  = sparse(nSo,nSo);              % iss <-> iss + ip <-> ip + ii <-> ii -- empty 
        case 3
            % only A{3} 
            A{1}  = sparse(nSo,nSo);
            A{2}  = sparse(nSo,nSo);
            A{3}  = ones(nSo,nSo)-eye(nSo);
        case 4
            % only A{3} with autoconnections                
            A{1}  = sparse(nSo,nSo);
            A{2}  = sparse(nSo,nSo);
            A{3}  = ones(nSo,nSo);
        case 5 
            % only B interactions
            A{1}  = sparse(nSo,nSo);
            A{2}  = sparse(nSo,nSo);
            A{3}  = sparse(nSo,nSo);
        case 6
            % autoconnections A{3}
            A{1}  = sparse(nSo,nSo);
            A{2}  = sparse(nSo,nSo);
            A{3}  = eye(nSo,nSo);
        case 7
            % autoconnections A + B
            A{1}  = eye(nSo,nSo);
            A{2}  = eye(nSo,nSo);
            A{3}  = eye(nSo,nSo);
        case 8
            % connections A
            A{1}  = ones(nSo,nSo);
            A{2}  = ones(nSo,nSo);
            A{3}  = ones(nSo,nSo);
        case 9
            % all autoconnections
            A{1}  = ones(nSo,nSo)-eye(nSo);
            A{2}  = ones(nSo,nSo)-eye(nSo);
            A{3}  = ones(nSo,nSo)-eye(nSo);
        case 10
            % autoconnections A + interactions B
            A{1}  = eye(nSo,nSo);
            A{2}  = eye(nSo,nSo);
            A{3}  = eye(nSo,nSo);
    
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
[pE,pC]     = spm_L_priors(nSo,pE,pC);      % spatial  priors
if ~isempty(customMode)
    [pE,pC]     = customMode(pE,pC,par);    % load custom priors
end
% Suppress channel noise
%--------------------------------------------------------------------------
pE.b        = pE.b - 16;
pE.c        = pE.c - 16;
% create LFP model
%--------------------------------------------------------------------------
M.dipfit.type = 'LFP';

M.IS        = 'spm_csd_mtf';
M.FS        = 'spm_fs_csd';
M.g         = 'spm_gx_erp';
%% alternative neural models
%% see -> spm_fx_cmc_tfm_gen    ( 8 neurons)  A{4} see slides - spm12/toolbox/NVC
%% see -> spm_fx_lfp            (13 neurons)
%% see -> spm_fx_fmri           ( 8 cells - BOLD)
%% see -> spm_fx_erp            ( 9 neurons)
M.f         = 'model_spm_fx_lfp';
M.x         = sparse(nSo,13);     
M.n         = nSo*13;
M.pE        = pE;
M.pC        = pC;
M.m         = nSo;            % areas
M.l         = nCh;            % electrodes        
M.Hz        = (1:64)';

M.Nmax      = 256;
M.nograph   = 1;              % no plot

% case nChannels > nSources
if nSources<nChannels
    M.pE.L=par.L;
    M.pC.L=par.L*mean(pC.L);
    M.g         = 'spm_gx_multiERP';
end


M_CSD       = M;