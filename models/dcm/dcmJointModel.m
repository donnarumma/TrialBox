function   DCM=dcmJointModel(data_trials,par)
% function DCM=dcmJointModel(data_trials,par)

% -> spm_csd_demo
%% Extract a data session
[~, Labels]         = getJointMonkeysLabels(1:24);
isel                = par.isdemo; %getSelectionIndexes(1,1:3);
InField             = par.InField;
%%
nTrials             = length(data_trials);
xY.y                = cell(1,nTrials);
xU.X                = nan(nTrials,length(Labels));
for iTrial=1:nTrials
    xY.y{iTrial}    = data_trials(iTrial).(InField);
    xU.X(iTrial,:)  = data_trials(iTrial).klabels;
end
NConds              = size(xU.X,2);         % number of conditions
nSources            = size(xY.y{1},2);      % number of Areas
nChannels           = nSources;             % number of total channels
%% test on subset of trials
iTrials             = any(xU.X(:,isel),2);
iConds              = sum(xU.X(iTrials,:))>0;
xU.X                = xU.X(iTrials,iConds);
xY.y                = xY.y(iTrials);
Labels              = Labels(iConds);
NConds              = size(xU.X,2);           % number of conditions
%%
xU.name             = Labels; 

%% MODEL
% M.IS = 'spm_csd_mtf'; % M.FS = 'spm_fs_csd'; % M.g  = 'spm_gx_erp'; % M.f  = 'spm_fx_lfp';
if ~isempty(par.custom_model)
    par.custom_model         =str2func(par.custom_model);
end
pms.nConditions = NConds;
pms.nSources    = nSources;
pms.nChannels   = nChannels;
pms.whichmodel  = par.whichmodel;
pms.customMode  = par.custom_model;
pms.donlfp      = par.donlfp;

M_CSD           = csdModel(pms);
sim_fun         = str2func(M_CSD.IS);
%%

%% invert - see: spm_dcm_csd
%--------------------------------------------------------------------------
[Ep,Cp,Eh,~,~,~,~,F_History] = DONNARUMMA_spm_nlsi_GN(M_CSD,xU,xY);

DCM.M        = M_CSD;             % Model ->
DCM.Hz       = M_CSD.Hz;          % frequency
%======================= RESULTS =======================
DCM.Ep       = Ep;                % conditional expectation
DCM.Cp       = Cp;                % conditional covariance
DCM.Ce       = exp(-Eh);          % ReML error covariance
% simulate learned spectral density 
%=======================       RESPONSE      ==============================
% [Hc, Hz] = spm_csd_mtf(Ep,M_CSD);
Hc           = sim_fun(Ep,M_CSD,xU);% [Hs, Hz, dtf] = spm_csd_mtf(Ep,M_CSD,DCM.xU);
DCM.Hc       = Hc;                  % conditional responses (y), channel space
Ec           = spm_unvec(spm_vec(xY.y) - spm_vec(Hc),Hc);     % prediction error
DCM.Rc       = Ec;                % conditional residuals (y), channel space
%==========================================================================
% predictions (source space - cf, a LFP from virtual electrode)
%--------------------------------------------------------------------------
if isfield(M_CSD,'U')
    M_Sou        = rmfield(M_CSD,'U');
else
    M_Sou        = M_CSD;
end
M_Sou.dipfit.type= 'LFP';

M_Sou.U      = 1; 
qp           = Ep;
qp.L         = ones(1,M_Sou.l);        % set virtual electrode gain to unity
qp.b         = qp.b - 32;              % and suppress non-specific and
qp.c         = qp.c - 32;              % specific channel noise

[Hs, ~, dtf] = spm_csd_mtf(qp,M_Sou,xU);
[ccf, pst]   = spm_csd2ccf(Hs,M_CSD.Hz);
[coh, fsd]   = spm_csd2coh(Hs,M_CSD.Hz);
DCM.dtf      = dtf;               % directed transfer functions (source space)
DCM.ccf      = ccf;               % cross  covariance functions (source space)
DCM.coh      = coh;               % cross  coherence  functions (source space)
DCM.fsd      = fsd;               % specific delay functions (source space)
DCM.pst      = pst;               % peristimulus time
DCM.Hs       = Hs;                % conditional responses (y), source space
%==========================================================================
%===================   spm_large_dcm_reduce   =============================
T       = full(spm_vec(M_CSD.pE));
sw      = warning('off','SPM:negativeVariance');
Pp      = spm_unvec(1 - spm_Ncdf(T,abs(spm_vec(Ep)),diag(Cp)),Ep);
Vp      = spm_unvec(full(diag(Cp)),Ep);
warning(sw);
DCM.Pp       = Pp;                % conditional probability
DCM.Vp       = Vp;                % conditional variances
%==========================================================================
                    %
DCM.F_History=F_History;

% DCM.A   = M_CSD.pE.A;
% DCM.B   = M_CSD.pE.B;
% DCM.C   = M_CSD.pE.C;

LFP         = cell(1,nTrials);
for iTrial=1:nTrials
    LFP{iTrial}=data_trials(iTrial).LFP';
end
DCM.xU      = xU;
DCM.xY.csd  = xY.y;
DCM.xY.lfp  = LFP;
DCM.Chamber = data_trials(1).Chamber;
DCM.xY.U    = [];
DCM.xY.Hz   = M_CSD.Hz;
DCM.xY.name={'monkey S'  'monkey K'};

return
nf          = ['LFP_M' num2str(whichmodel) '_' session_name add];
fn =[save_dir nf];
fprintf('Saving %s\n',fn);
DCM.nf=nf;
DCM.fn=fn;
DCM.par=par;
save(fn,'DCM');


for iTrial=1:length(DCM.xY.csd)
    fprintf('Trial %g\n',iTrial);
    trialError(iTrial)=computeError(DCM.xY.csd{iTrial},Hc{iTrial});
    disp(trialError(iTrial));
end
DCM.trialError = trialError;

if ifserver 
    return;
else
    DCM_RESULTS_LFP_KS(DCM)
    %tentativePlots(DCM);
end
return

params.hfig = [];
hCSD=DONNARUMMA_plotCSD(M_CSD.Hz,DCM.Hc{iTrial},DCM.xY.y{iTrial},params);

params.hfig=[];
hF=DONNARUMMA_plotFreeEnergy(DCM.F_History,params);   

params.hfig=[]; P_TRUE=[];
hEp=plot_conditionalExpectation(DCM.Ep,P_TRUE,params);

CM.A=DCM.Ep.A;
CM.C=DCM.Ep.C;
if ~isempty
    AT.A=P_TRUE.A;
    AT.C=P_TRUE.C;
else
    AT=[];
end
params.hfig=[];
hEpA=plot_conditionalExpectation(CM,AT,params);


%%
plot_spm_dcm_csd_results(DCM,'Coupling (A)',figure);
plot_spm_dcm_csd_results(DCM,'spectral data',figure);  
plot_spm_dcm_csd_results(DCM,'Coupling (B)',figure);
plot_spm_dcm_csd_results(DCM,'Coupling (C)',figure);
plot_spm_dcm_csd_results(DCM,'trial-specific effects',figure);
plot_spm_dcm_csd_results(DCM,'Input',figure);
plot_spm_dcm_csd_results(DCM,'Transfer functions',figure);
plot_spm_dcm_csd_results(DCM,'Cross-spectra (sources)',figure)
plot_spm_dcm_csd_results(DCM,'Cross-spectra (channels)',figure)
plot_spm_dcm_csd_results(DCM,'Coherence (sources)',figure)
plot_spm_dcm_csd_results(DCM,'Coherence (channels)',figure)
plot_spm_dcm_csd_results(DCM,'Covariance (sources)',figure)
plot_spm_dcm_csd_results(DCM,'Covariance (channels)',figure)
% plot_spm_dcm_csd_results(DCM,'Dipoles',figure);
%%
