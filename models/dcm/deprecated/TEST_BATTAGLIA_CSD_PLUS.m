%% function TEST_BATTAGLIA_CSD_PLUS
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
session_name                    = 'SK004';%{'SK001','SK009'};  % session name 
idir                            = 1;        % directions -> 1-8 
% fixed parameteres
nLFP                            = 5;        % max number of LFPs
par.BATTAGLIA_CSD               = BATTAGLIA_CSDParams();
S                               = filesep;
if isonserver
    par.BATTAGLIA_CSD.save_dir      = ['~' S 'SAPIENZA' S 'DCM_SERVER' S session_name S 'dir' num2str(idir) S];
else 
    par.BATTAGLIA_CSD.save_dir      = ['~' S 'TESTS' S 'SAPIENZA' S 'DCM_LOCAL' S session_name S 'dir' num2str(idir) S];
end
if ~isfolder(par.BATTAGLIA_CSD.save_dir)
    mkdir(par.BATTAGLIA_CSD.save_dir);
end
fprintf('Saving in %s\n',par.BATTAGLIA_CSD.save_dir);
par.BATTAGLIA_CSD.session_name  = session_name;
par.BATTAGLIA_CSD.whichmodel    = 7;        % Self Model   


session_data = BattagliaArrangeTrials(par);
%% Monkey S Action model prior
fprintf('Monkey S Action model prior\n')
par.BATTAGLIA_CSD.isdemo        = getSelectionIndexes(idir,1);  % get S Action Trials (1)
par.BATTAGLIA_CSD.selK          = [];                           % do not take K LFP
DCM_S                           = cell(0,0);
for iLFP=1:nLFP
    %fprintf('Trying iLFP %g\n',iLFP)
    par.BATTAGLIA_CSD.selS      = iLFP;
    DCM_cur                     = BATTAGLIA_CSD(par.BATTAGLIA_CSD);
    if ~isempty(DCM_cur)
        DCM_S{iLFP}             = DCM_cur;
    else
        fprintf('Maximum iLFP reached %g\n',iLFP-1)
        break
    end
end
%pause
%% Monkey K Action model prior
fprintf('Monkey K Action model prior\n')
par.BATTAGLIA_CSD.isdemo        = getSelectionIndexes(idir,2); % get K Action Trial (2)
par.BATTAGLIA_CSD.selS          = [];   % do not take S LFP
DCM_K                           = cell(0,0);
for iLFP=1:nLFP
    %fprintf('Trying iLFP %g\n',iLFP)
    par.BATTAGLIA_CSD.selK      = iLFP;
    DCM_cur                     = BATTAGLIA_CSD(par.BATTAGLIA_CSD);
    if ~isempty(DCM_cur)
        DCM_K{iLFP}             = DCM_cur;
    else
        fprintf('Maximum iLFP reached %g\n',iLFP-1);
        break;
    end
end
%pause
%% select winner Priors
for iS=1:length(DCM_S)
    Free_S(iS)=max(DCM_S{iS}.F_History);
end
[~,S_winner]=max(Free_S);

for iK=1:length(DCM_K)
    Free_K(iK)=max(DCM_K{iK}.F_History);
end
[~,K_winner]=max(Free_K);
fprintf('S LFP index %g\n', S_winner);
fprintf('K LFP index %g\n', K_winner);
%% Joint Model
%% F=1.683e+04 err=7.193e-03
%% ~/TESTS/SAPIENZA/DCM_SERVER/LFP_M5_SK009_D7_S2_K1_C_1_9_17_customPriorsv3
fprintf('Monkey S-K model\n')
DCM_J                         = cell(0,0);
par.BATTAGLIA_CSD.custom_model= 'customPriorsv3';
par.BATTAGLIA_CSD.donlfp      = false;
par.BATTAGLIA_CSD.isdemo      = getSelectionIndexes(idir,1:3); % get S, K and S-K Trials
par.BATTAGLIA_CSD.S_PRIORS    = DCM_S{S_winner}.fn;
par.BATTAGLIA_CSD.K_PRIORS    = DCM_K{K_winner}.fn;
par.BATTAGLIA_CSD.ifserver    = true; %false;  % do not plot result on server
par.BATTAGLIA_CSD.whichmodel  = 5;      % Interaction Model 
par.BATTAGLIA_CSD.selS        = S_winner;
par.BATTAGLIA_CSD.selK        = K_winner;
DCM_J{end+1}                      = BATTAGLIA_CSD(par.BATTAGLIA_CSD);
%% ~/TESTS/SAPIENZA/DCM_SERVER/LFP_M5_SK009_D7_S2_K1_C_1_9_17_customPriorsv3_DONLFP
%par.BATTAGLIA_CSD.donlfp      = true;
%DCM_J{end+1}                      = BATTAGLIA_CSD(par.BATTAGLIA_CSD);
%% ~/TESTS/SAPIENZA/DCM_SERVER/SK001/dir1/LFP_M5_SK001_D7_S2_K1_C_1_9_17
par.BATTAGLIA_CSD.donlfp      = false;
par.BATTAGLIA_CSD.custom_model= [];
DCM_J{end+1}                  = BATTAGLIA_CSD(par.BATTAGLIA_CSD);
%% ~/TESTS/SAPIENZA/DCM_SERVER/LFP_M5_SK009_D7_S2_K1_C_1_9_17_DONLFP
%par.BATTAGLIA_CSD.donlfp      = true;
%par.BATTAGLIA_CSD.custom_model= [];
%DCM_J{end+1}                  = BATTAGLIA_CSD(par.BATTAGLIA_CSD);
%% Winner J err=1.456e-04 FreeEnergy=1.661e+04
for iJ=1:length(DCM_J)
    Free_J(iJ)=max(DCM_J{iJ}.F_History);
end
[~,J_winner]=max(Free_J);
fprintf('Joint model winner: %g\n',J_winner)
fprintf('Chamber %s\n', DCM_J{J_winner}.Chamber);
DCM_RESULTS_LFP_KS(DCM_J{J_winner})