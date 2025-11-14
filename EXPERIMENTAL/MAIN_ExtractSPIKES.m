%% MAIN_ExtractSPIKES
% copia della funzione S_translation di Donnarumma

%clear; close all
% cd '.......'
filesep = '/'; 

%addpath(genpath('/home/frosolone/EXPERIMENT_DCM_RESULTS/Extract_Spikes/OnedriveFra'));
addpath("GOODCELLS\")
addpath("SKH_SUA_frontal_192\")
interval = [200,1000]; % Time max overall 3.1 s or 2.2

selMonkey = 'all'; %%% either S, K (?) or all


% All Frontal Sessions...26 o 27?? la 26 non ci sta...la 27 sì
% session_list = {'SK011','SK012','SK021','SK022','SK023','SK024','SK025','SK027','SK028','SK030','SK032',...};
%                 'SK033','SK035','SK036','SK038','SK039','SK040','SK041','SK042','SK043','SK044','SK045','SK046','SK047','SK048'};

% Frontal Sessions with 192 Trials
session_list = {'SK043'}
  %  'SK022','SK033','SK025','SK035','SK036','SK038','SK042'};

VARIABLES_monkey_joint_action;
try
    MODE_ALIGN;
catch
    %     MODE_ALIGN=MOVEMENT_ONSET_MODE;
    MODE_ALIGN=TARGET_APPEARS_MODE;
end

S_sessions = struct();
for n_sess = 1:length(session_list)
    session_name = session_list{n_sess};
    try
        session_name = strcat(session_name,'H_SUA');
    catch
        session_name='SK033H_SUA.mat';
        fprintf('S_translation Warning: using %s\n',session_name);%pause;
    end
    d=load (session_name);
    Trials=d.Trials;
    chamber = Trials(1).Chamber;
    ti=tic;
    fprintf('Getting data from %s\n',session_name);
    S=    getSpikesJointAllNeurons(Trials);
    fprintf('Elapsed time: %g s\n',toc(ti))
    S=      getSuccessfullTrials(S,Trials);
    S=   getJointConditionTrials(S,Trials);
    S= getSolo_S_ConditionTrials(S,Trials);
    S= getSolo_K_ConditionTrials(S,Trials);
    S=             getAreaTrials(S,Trials);
    S=            getGoodNeurons(S,Trials);
    % S=      getTuningSoloNeurons(S,Trials);
    % S=     getTuningJointNeurons(S,Trials);
    % S=       getTuningObsNeurons(S,Trials);
    % S=    getTuningNoSoloNeurons(S,Trials);
    % S=   getTuningNoJointNeurons(S,Trials);
    % S=     getTuningNoObsNeurons(S,Trials);
    % if length(CELL_SELECTION)>1
    %     S=          getTuningNeurons(S,Trials,CELL_SELECTION(2));
    %     S=        getNoTuningNeurons(S,Trials,CELL_SELECTION(2));
    % end
    S=getTimeMovementOnsetTargetTrials(S,Trials); %% new entry
    for ind=1:8
        S=getTConditionTrials(S,Trials,ind);
    end
    S=       getTimeTargetTrials(S,Trials);
    %S=retuneSpikeWithTimeStamp1(S,Trials);
    if MODE_ALIGN==TARGET_APPEARS_MODE
        S=retuneSpikeWithTargetAppears(S,Trials);
    elseif MODE_ALIGN==MOVEMENT_ONSET_MODE
        S=retuneSpikeWithMovementOnset(S,Trials);
    elseif MODE_ALIGN==MOVEMENT_ONSET_K_MODE
        S=retuneSpikeWithMovementOnsetK(S,Trials);
    elseif MODE_ALIGN==MOVEMENT_ONSET_S_MODE
        S=retuneSpikeWithMovementOnsetS(S,Trials);
    end
    % solo_S=S.solo_S_trials.*S.successfull_trials;
    % solo_K=S.solo_K_trials.*S.successfull_trials;
    % joint = S.joint_trials.*S.successfull_trials;

    try
        interval;
    catch
        interval=[-500,1000];
        interval=[-100,300];

        interval = [-199,400];
    end

    try
        selMonkey;
    catch
        selMonkey='all';
    end
    S=prepareSpikeAndTrials(S,Trials,interval,selMonkey);    %v1
    % S=prepareSpikeAndTrialsPlus(S,Trials,interval,selMonkey,CELL_SELECTION(1));
    % %v2 % non funziona

    dat=S.dat;

    % t_S_on = [S.time_S_movement_onset];
    % t_S_on = t_S_on(logical([S.successfull_trials]));
    %
    % t_K_on = [S.time_K_movement_onset];
    % t_K_on = t_K_on(logical([S.successfull_trials]));
    %
    % t_targApp = [S.time_target_trials];
    % t_targApp = t_targApp(logical([S.successfull_trials]));
    %
    % for iTr=1:length(dat)
    %     dat(iTr).time_S_Movement_Onset = t_S_on(iTr);
    %     dat(iTr).time_K_Movement_Onset = t_K_on(iTr);
    %     dat(iTr).time_Target_Appears = t_targApp(iTr);
    % end

    S_sessions(n_sess).data = dat;
    for iTr=1:length(dat)
        S_sessions(n_sess).data(iTr).name = session_list{n_sess};
    end
end

t = linspace(interval(1),interval(2),size(dat(1).spikes,2));

data_session = [S_sessions(:).data];

data_trials = struct();

for iTr=1:length(data_session)
    %iTr
    data_trials(iTr).trialId                = iTr;
    data_trials(iTr).timeSpikes             = t;
    data_trials(iTr).Spikes                 = data_session(iTr).spikes;
    data_trials(iTr).trialTypeDir           = data_session(iTr).trialType;
    data_trials(iTr).trialTypeCond          = data_session(iTr).trialType2;
    data_trials(iTr).OldTrialId             = data_session(iTr).trialId;
    data_trials(iTr).SessName               = data_session(iTr).name; 
    data_trials(iTr).chamber                = chamber;


    % data_trials(iTr).time_S_Movement_Onset  = dat(iTr).time_S_Movement_Onset;
    % data_trials(iTr).time_K_Movement_Onset  = dat(iTr).time_K_Movement_Onset;
    % data_trials(iTr).time_Target_Appears    = dat(iTr).time_Target_Appears;
end

data_trials(1).All_N_Identity         = S.AllNeuronIdentity
data_trials(1).All_Good_N_Identity    = S.AllGoodNeuronIdentity
%% SAVE DATA_TRIALS - OPTIONAL
save_path = pwd;
save_path = strcat(save_path,filesep,'SpikesExtracted',filesep);

if ~exist(save_path, 'dir')
    mkdir(save_path);
end

name_save = strcat(selMonkey,'_',num2str(interval(1)),'_',num2str(interval(2)),'_',S.area, session_name);
save(fullfile(save_path, strcat(name_save,'.mat')), 'data_trials');