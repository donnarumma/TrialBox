%% function TEST_BATTAGLIA_DCM
clear all
rng(10)
%% check if is on server
[~,nn]=system('hostname'); nn=strtrim(nn);
if strcmp(nn,'rk018610')
    isonserver=true;
else
    isonserver=false;
end
if ~isonserver
    test_dir    ='~/TESTS/SAPIENZA/DCM';
else
    test_dir    ='~/SAPIENZA/SERVER/DCM';
end

% fixed parameters
session_name                    = 'SK009';%'SK004';%{'SK001','SK009'};  % session name 
idir                            = 1;        % directions -> 1-8 
S                               = filesep;
%% Step 0: arrange trials
fprintf('Step 0: arrange trials\n');
par.BattagliaArrangeTrials              = BattagliaArrangeTrialsParams();
% par.BattagliaArrangeTrials.whichmodel   = 7;        % 7 Model for action session. 5 Model for all sessions
par.BattagliaArrangeTrials.isdemo       = 1;        % getSelectionIndexes(1,1:3);
par.BattagliaArrangeTrials.selS         = 2;
par.BattagliaArrangeTrials.selK         = 1;
par.BattagliaArrangeTrials.session_name = session_name;  % which session
data_trials                             = BattagliaArrangeTrials(par.BattagliaArrangeTrials);
%% Step 1: CSD compute
fprintf('Step 1: csd compute\n');
signal_process                          = 'CSD';
par.csdCompute                          = csdComputeParams();            
par.csdCompute.InField                  = 'LFP';
par.csdCompute.OutField                 = signal_process;
par.csdCompute.rsfactor                 = 0.2;
par.csdCompute.optrescale               = 3;%-1;
data_trials                             = csdCompute(data_trials,par.csdCompute);
test_dir                                = [test_dir '_opt' num2str(par.csdCompute.optrescale) S];
%% learning
fprintf('Step 2: Monkey S-K model\n')
par.dcmJointModel.whichmodel            = 5;
par.dcmJointModel.custom_model          = [];%'customPriorsv3';
par.dcmJointModel.InField               = signal_process;
par.dcmJointModel.donlfp                = false;
par.dcmJointModel.isdemo                = getSelectionIndexes(idir,1:3); % get S, K and S-K Trials
par.dcmJointModel.S_PRIORS              = [];%DCM_S{S_winner}.fn;
par.dcmJointModel.K_PRIORS              = [];%DCM_K{K_winner}.fn;
out.dcmJointModel.DCM                   = dcmJointModel(data_trials,par.dcmJointModel);
%% save test 
% custom directory and file name string construction
[~,Labels]  = getJointMonkeysLabels(1:24);
iConds       = find(ismember(Labels,DCM_Joint.xU.name));
save_dir    =[test_dir session_name S];
save_dir    =[save_dir 'Dir'];
for ind=1:length(idir)
    save_dir=[save_dir '_' num2str(idir(ind))];
end
save_dir    =[save_dir S];
save_file   =['LFP_M' num2str(par.dcmJointModel.whichmodel) '_' session_name];
save_file   =[save_file '_D' num2str(par.BattagliaArrangeTrials.dmode)];
save_file   =[save_file '_S' num2str(par.BattagliaArrangeTrials.selS)];
save_file   =[save_file '_K' num2str(par.BattagliaArrangeTrials.selK)];

save_file=[save_file '_C'];
for inc=1:length(iConds)
    save_file=[save_file '_' num2str(iConds(inc))];
end
if ~isempty (par.dcmJointModel.custom_model)
    save_file=[save_file '_' par.dcmJointModel.custom_model];
end
if par.dcmJointModel.donlfp
    save_file = [save_file '_DONLFP'];
end
fprintf('Saving in %s\n',[save_dir save_file]);
DCM         =out.dcmJointModel.DCM;
DCM.file    =save_file;
DCM.fullpath=[save_dir, save_file];
DCM.par     =par.dcmJointModel;
if ~isfolder(save_dir)
    mkdir(save_dir);
end
save(DCM.fullpath,'DCM','data_trials');
%% stop if is on server
if isonserver; return; end
%% LFP
iTrial          = 1;
fprintf('Showing Trial %g\n',iTrial);
hfg.lfp         = figure;
hold on; box on; grid on;
lw              = 2;
InField         = 'LFP';
xfld            = 't';
TField          = [xfld InField];
col1            = [0,0,0];
nSources    = size(data_trials(iTrial).(InField),1);
tiledlayout(nSources,1);
sigtime     = data_trials(iTrial).(TField);
sigtime     = sigtime-sigtime(1);
for iSource=1:nSources
    nexttile;
    hold on; box on; grid on;
    plot(sigtime,data_trials(iTrial).(InField)(iSource,:),'color',col1,'linewidth',lw);
    xlabel('time [s]')
    ylabel('LFP [mV]') 
    title(['Source ' num2str(iSource)]);
end
sgtitle(['Low Field Potential, Trial ' num2str(iTrial)])
%% CSD
iTrial          = 2;
fprintf('Showing CSD result of Trial %g\n',iTrial);
hfg.cross       = figure;
col1            = [1,0,0];
col2            = [0,1,0];
lw              = 3;
InField         = 'CSD';
xfld            = 't';
TField          = [xfld,InField];
[~,nRows,nCols] = size(data_trials(iTrial).(InField));
p               = [0,0,1500,600];
fs              = 13;
tiledlayout(nRows,nCols);
for iRow=1:nRows
    for iCol=1:nCols
        nexttile;
        hold on; box on; grid on;
        % plot original
        xl  = [data_trials(iTrial).(TField)(1),data_trials(iTrial).(TField)(end)];
        plot(data_trials(iTrial).(TField), ...
             real(data_trials(iTrial).(InField)(:,iRow,iCol)),'color',col1,'linewidth',lw,'linestyle',':');%,M.Hz,real(CSD(:,1,2)),':')
        % plot reconstructed
        plot(data_trials(iTrial).(TField), ...
             real(DCM.Hc{iTrial}(:,iRow,iCol)),'color',col2,'linewidth',lw);%,M.Hz,real(CSD(:,1,2)),':')
        xlim(xl)
        title(['Source ' num2str(iRow) ' vs Source' num2str(iCol)])
        if iCol==nCols && iRow==nRows
            legend('original','reconstructed');
        end
        if iRow==nRows
            xlabel('Freq [Hz]');
        end
        if iCol==1
            ylabel('CSD [dB/Hz]');
        end
        set(gca,'fontsize',fs);
    end
end
trialError  = computeError(data_trials(iTrial).(InField),DCM.Hc{iTrial});
st          = sprintf('[cross]-spectral density - RMSE: %0.3f', trialError.RMSE);
sgtitle(st,'fontsize',fs+2)
set(hfg.cross,'Position',p);
set(hfg.cross, 'PaperPositionMode','auto');

%%
% plot parameters and estimates
%--------------------------------------------------------------------------
hfg.expects = figure;
hold on; box on; grid on;
bar(exp(spm_vec(DCM.Ep)))
title('conditional expectation')
errorbar(exp(spm_vec(DCM.Ep)),exp(DCM.Cp))

%% some plots
plot_spm_dcm_csd_results(DCM,'Coupling (B)',figure);
plot_spm_dcm_csd_results(DCM,'trial-specific effects',figure);

% DCM_RESULTS_LFP_KS(DCM);
%%
