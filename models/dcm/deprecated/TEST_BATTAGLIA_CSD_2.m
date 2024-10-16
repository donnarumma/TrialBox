%% function TEST_BATTAGLIA_CSD_2
clear all
rng(10)
%% check if is on server
[~,nn]=system('hostname'); nn=strtrim(nn);
if strcmp(nn,'rk018610')
    isonserver=true;
else
    isonserver=false;
end
%%
session_name                    = 'SK009';%SK004';%{'SK001','SK009'};  % session name 
idir                            = 1;        % directions -> 1-8 
% fixed parameteres
nLFP                            = 5;        % max number of LFPs
par.BATTAGLIA_CSD               = BATTAGLIA_CSDParams();
par.BATTAGLIA_CSD.save_dir      = '';
par.BATTAGLIA_CSD.session_name  = session_name;
par.BATTAGLIA_CSD.whichmodel    = 7;        % Self Model   
%% Joint Model
%% F=1.683e+04 err=7.193e-03
%% ~/TESTS/SAPIENZA/DCM_SERVER/SK001/dir1/LFP_M5_SK001_D7_S2_K1_C_1_9_17
par.BATTAGLIA_CSD.custom_model= [];
par.BATTAGLIA_CSD.donlfp      = false;
par.BATTAGLIA_CSD.isdemo      = getSelectionIndexes(idir,1:3); % get S, K and S-K Trials
par.BATTAGLIA_CSD.S_PRIORS    = [];
par.BATTAGLIA_CSD.K_PRIORS    = [];
par.BATTAGLIA_CSD.ifserver    = true; %false;  % do not plot result on server
par.BATTAGLIA_CSD.whichmodel  = 5;      % Interaction Model 
par.BATTAGLIA_CSD.selS        = 2;
par.BATTAGLIA_CSD.selK        = 1;
par.BATTAGLIA_CSD.donlfp      = false;
DCM_Joint                     = BATTAGLIA_CSD(par.BATTAGLIA_CSD);
%% ~/TESTS/SAPIENZA/DCM_SERVER/SK001/dir1/LFP_M5_SK009_D7_S2_K1_C_1_9_17_DONLFP
DCM_RESULTS_LFP_KS(DCM_Joint);
return
%% Winner J err=1.456e-04 FreeEnergy=1.661e+04
for iJ=1:length(DCM_J)
    Free_J(iJ)=max(DCM_J{iJ}.F_History);
end
[~,J_winner]=max(Free_J);
fprintf('Joint model winner: %g\n',J_winner)
fprintf('Chamber %s\n', DCM_J{J_winner}.Chamber);
DCM_RESULTS_LFP_KS(DCM_J{J_winner})