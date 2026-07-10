%% MAIN_BEHAVIOUR_SPIKES_EXTRACT
% MAIN_BEHAVIOUR_SPIKES_EXTRACT.m
%
% Extracts behavioural metrics from SUA sessions and aggregates results
% across all selected sessions.
%
% Main outputs:
%   - RT
%   - ST
%   - PV
%   - PVT
%   - Trajectory metrics
%
% Output:
%   Data.Behav structure saved as MAT file.
%
% Author: Mirco Frosolone
% Date: 2026-07-06
%

clear; close all;
par.irng = 10;
rng(par.irng)

params = struct();

% default params
defParams = getDefaultSPIKESparams;

% merge params
params = mergeParams(defParams, params);

%% Local Params
warning('Verify the extraction interval before running the script.');

interval = [0, 0.5];


% FRONTAL
session_list        = {...
    'SK007','SK008','SK009','SK010','SK011','SK012','SK013','SK014', ...
    'SK020','SK021','SK022','SK023','SK024','SK025','SK027','SK028', ...
    'SK029','SK030','SK031','SK032','SK033','SK034','SK035','SK036', ...
    'SK037','SK038','SK039','SK040','SK041','SK042','SK043','SK044', ...
    'SK045','SK046','SK047','SK048'};
params.chamber = 'Frontal';

% 
% %% PARIETAL
% session_list = {...
%     'SK003','SK004','SK049','SK050','SK051','SK052','SK053','SK054', ...
%     'SK055','SK056','SK057','SK058','SK059','SK060','SK061','SK062', ...
%     'SK063','SK064','SK065','SK066','SK067','SK069','SK070','SK071', ...
%     'SK072' };
% params.chamber = 'Parietal';


directory_save  = 'BEHAVIOUR_DATA';
file_name       = ['Behaviour_',params.chamber,'_TotalSessionsData_',...
    num2str(interval(1)),'_',erase(num2str(interval(2)),"."),'.mat'];

fprintf('\n');
fprintf('================ BEHAVIOUR SPIKES EXTRACTION PARAMETERS ================\n');
fprintf('Chamber        : %s\n', params.chamber);
fprintf('Time Interval  : [%g  %g] ms\n', ...
        interval(1)*1000,...
        interval(2)*1000);
% fprintf('Mode Aligned   : %s\n', modeAlignToString(params.MODE_ALIGN));
fprintf('Monkey         : %s\n', params.selMonkey);
fprintf('N. Sessions    : %d\n', numel(session_list));
fprintf('===============================================================\n\n');

%% Step 0: arrange trials of all sessions
% Note: selS and selK are related to LFP analyses and are not used here.
RT_session          = struct();
PV_session          = struct();
PVT_session         = struct();
ST_session          = struct();
Jxy_session         = struct();

for idsess = 1:numel(session_list)

    session_name = [session_list{idsess}, 'H_SUA'];

    try
        d = load(session_name);
    catch ME
        warning('Unable to load session %s: %s', ...
            session_name, ME.message);
        continue
    end

    Trials  = d.Trials;
    [~, ~, klabels, success] = getClassInfo(Trials);
    data_SUA     = Trials(success);
    klabels      = klabels(success,:);

    tstartExtract       = -0.7;
    tendExtract         =  5.0;
    AlignEvent          = 2024;% Align to peripheral target onset
    LockEvent           = 2024;
    for i=1:numel(data_SUA)
        TimeLock=data_SUA(i).Event (data_SUA(i).Event (:,1)==AlignEvent,2); %Get relative time of event presentation
        NewTime=linspace(data_SUA(i).TimeStamp(1),data_SUA(i).TimeStamp(2),size(data_SUA(i).JXYEXY,2));
        [StartP,~] =Odysseas_Get_Nearest_Value(tstartExtract+TimeLock ,NewTime);%Get start centered on the event
        [StopP,~] =Odysseas_Get_Nearest_Value(tendExtract+TimeLock,NewTime);  %Get stop centered on the event
        data_SUA(i).JXYEXY=data_SUA(i).JXYEXY(:,StartP:StopP); %Get the newly centered matrix
        data_SUA(i).timeJXYEXY = NewTime(StartP:StopP);
        tmp = find(klabels(i,:),1,'first');
        data_SUA(i).trialType = tmp;
        data_SUA(i).klabels = klabels(i,:);
    end
    data_SUA       = RTextract(data_SUA);
    data_SUA       = ETextract(data_SUA);
    

    par.alignBattaglia.AlignEvent   = AlignEvent;
    par.alignBattaglia.LockEvent = LockEvent;
    par.alignBattaglia.InField      = 'JXYEXY';
    par.alignBattaglia.OutField     = 'JXYEXY';

    data_SUA                     = alignBattaglia(data_SUA,par.alignBattaglia);

    dataInField                     = 'JXYEXY';
    dataxfld                        = 'time';
    tStart                  = interval(1);
    tEnd                    = interval(2);
    for iTrial=1:numel(data_SUA)
        time_app            = data_SUA(iTrial).([dataxfld dataInField]);
        zero_ind            = find(time_app >= 0, 1, 'first');

        if isempty(zero_ind)
            error('No non-negative sample found in time vector.');
        end

        tStart_ind          = find(time_app >= (time_app(zero_ind) + tStart), 1, 'first');
        tEnd_ind            = find(time_app >= (time_app(zero_ind) + tEnd), 1, 'first');
        if isempty(tEnd_ind)
            tEnd_ind = length(time_app);
        end
        data_SUA(iTrial).(dataInField)               = data_SUA(iTrial).(dataInField)(:, tStart_ind:tEnd_ind);
        data_SUA(iTrial).([dataxfld dataInField])    = data_SUA(iTrial).([dataxfld dataInField])(tStart_ind:tEnd_ind);
    end

    [~,Labels] = getJointMonkeysLabels(1:24);

    shiftRawSUA = zeros(8, 1);
    shiftRawSUA(1:2) = [0.5000; -3.8000]; % Offset per S
    shiftRawSUA(5:6) = [1.0000; -3.0000]; % Offset per K

    for iTrial=1:length(data_SUA)
        data_SUA(iTrial).trialId = iTrial;
        data_SUA(iTrial).trialName        = Labels{data_SUA(iTrial).trialType};
        data_SUA(iTrial).JXYEXY = data_SUA(iTrial).JXYEXY - shiftRawSUA;
    end
    %% Step 1: Arrange data
    data_trials = findSKdirection(data_SUA);

    num_dir         = 8; % number of direction
    num_cond        = 3; % number of conditions: 1-SoloS 2-SoloK 3-Joint S-K


    Condition = cell(num_cond,num_dir);
    for cd=1:num_cond
        for ndir=1:num_dir
            Condition{cd,ndir} = find([data_trials.Condition]==cd & [data_trials.Direction]==ndir);
        end
    end

    Move = struct();
    lenM = NaN(num_cond,num_dir);
    for cd=1:num_cond
        for ndir=1:num_dir
            Move(cd).(strcat('dir',num2str(ndir)))  = data_trials(Condition{cd,ndir});
            lenM(cd,ndir)                           = length(data_trials(Condition{cd,ndir}));
        end
    end

    %% Trajectory Error
    par.TrajXYextrat4eval.num_cond          = num_cond;
    par.TrajXYextrat4eval.num_dir           = num_dir;
    par.TrajXYextrat4eval.fieldname         = 'JXYEXY';
    par.TrajXYextrat4eval.time              = 'timeJXYEXY';
    [JxyS,JxyK,Jxy_name,Jxy_time,RT_Jxy]    = TrajXYextrat4eval(Move,par.TrajXYextrat4eval);
    [Jxy_trajS,Jxy_trajK]                   = XYtrajectoryEval(JxyS,JxyK);
    Jxy_session(idsess).session             = session_list{idsess};
    Jxy_session(idsess).S                   = JxyS;
    Jxy_session(idsess).K                   = JxyK;
    Jxy_session(idsess).name                = Jxy_name;
    Jxy_session(idsess).time                = Jxy_time;
    Jxy_session(idsess).RT_Jxy              = RT_Jxy;
    Jxy_session(idsess).TrajS               = Jxy_trajS;
    Jxy_session(idsess).TrajK               = Jxy_trajK;


    par.RTextrat4eval.num_cond  = num_cond;
    par.RTextrat4eval.num_dir   = num_dir;
    par.RTextrat4eval.idsess    = idsess;
    par.RTextrat4eval.fieldname = 'RT';

    [RT_tot,~,~,~,~]       = RTextrat4eval(Move,lenM,par.RTextrat4eval);
    RT_session(idsess).session  = session_list{idsess};
    RT_session(idsess).cond     = RT_tot;

    par.PVExtrat4eval.num_cond       = num_cond;
    par.PVExtrat4eval.num_dir        = num_dir;
    par.PVExtrat4eval.idsess         = idsess;
    par.PVExtrat4eval.fieldname      = 'JXYEXY';
    par.PVExtrat4eval.time           = 'timeJXYEXY';
    par.PVExtrat4eval.windowSize     = 50;

    [PV_tot,~,~,~,~]                = PVExtrat4eval(Move,lenM,par.PVExtrat4eval);
    PV_session(idsess).session              = session_list{idsess};
    PV_session(idsess).cond                 = PV_tot;

    par.PVTExtrat4eval.num_cond            = num_cond;
    par.PVTExtrat4eval.num_dir             = num_dir;
    par.PVTExtrat4eval.idsess              = idsess;
    par.PVTExtrat4eval.fieldname           = 'JXYEXY';
    par.PVTExtrat4eval.time                = 'timeJXYEXY';
    par.PVTExtrat4eval.windowSize          = 50;


    [PVT_Tot,~,~,~,~]              = PVTExtrat4eval(Move,lenM,par.PVTExtrat4eval);
    PVT_session(idsess).session             = session_list{idsess};
    PVT_session(idsess).cond                = PVT_Tot;

    par.STextrat4eval.num_cond              = num_cond;
    par.STextrat4eval.num_dir               = num_dir;
    par.STextrat4eval.idsess                = idsess;
    par.STextrat4eval.fieldname             = 'ET';

    [ST_tot,~,~,~,~]                   = STextrat4eval(Move,lenM,par.STextrat4eval);
    ST_session(idsess).session              = session_list{idsess};
    ST_session(idsess).cond                 = ST_tot;
end
%% Trajectory parameters
par.computeTrajectoryExitDir.r_dim = 1.81;
Jxy_Result = computeTrajectoryExitDir(Jxy_session,par.computeTrajectoryExitDir);

[Jxy_S, Jxy_K] = JxyMatrix(Jxy_Result);

par.BehavWinner_Trajectory.num_dir              = num_dir;
par.BehavWinner_Trajectory.InField              = 'angular_error';
[Jxy_WinnerAE,Jxy_TabSoloAE,Jxy_TabJointAE,Jxy_compareAE]  = BehavWinner_Trajectory(Jxy_S,Jxy_K, par.BehavWinner_Trajectory);
par.BehavWinner_Trajectory.InField              = 'angular_error_magnitude';
[Jxy_WinnerAE_abs,Jxy_TabSoloAE_abs,Jxy_TabJointAE_abs,Jxy_compareAE_abs]  = BehavWinner_Trajectory(Jxy_S,Jxy_K, par.BehavWinner_Trajectory);
par.BehavWinner_Trajectory.InField              = 'exit_time';
[Jxy_WinnerET,Jxy_TabSoloExT,Jxy_TabJointExT,Jxy_compareExT]       = BehavWinner_Trajectory(Jxy_S,Jxy_K, par.BehavWinner_Trajectory);
par.BehavWinner_Trajectory.InField              = 'ExitCorr';
[Jxy_WinnerEC,Jxy_TabSoloEC,Jxy_TabJointEC,Jxy_compareEC]          = BehavWinner_Trajectory(Jxy_S,Jxy_K, par.BehavWinner_Trajectory);
par.BehavWinner_Trajectory.InField              = 'CMT';
[Jxy_WinnerCMT,Jxy_TabSoloCMT,Jxy_TabJointCMT,Jxy_compareCMT]      = BehavWinner_Trajectory(Jxy_S,Jxy_K, par.BehavWinner_Trajectory);

par.behavSpatialWinner.num_dir              = num_dir;
Jxy_Table_Spatial = behavSpatialWinner(Jxy_compareAE,par.behavSpatialWinner);
%% MEAN FOR SESSION AND DIRECTION Trajectory
par.BehavWinner_TrajectoryMEAN.num_dir              = num_dir;
par.BehavWinner_TrajectoryMEAN.InField              = 'angular_error';
[Jxy_WinnerMeanAE,Jxy_TabSoloMeanAE,Jxy_TabJointMeanAE]          = BehavWinner_Trajectory_MEAN(Jxy_S,Jxy_K, par.BehavWinner_TrajectoryMEAN);
par.BehavWinner_TrajectoryMEAN.InField              = 'angular_error_magnitude';
[Jxy_WinnerMeanAE_abs,Jxy_TabSoloMeanAE_abs,Jxy_TabJointMeanAE_abs]          = BehavWinner_Trajectory_MEAN(Jxy_S,Jxy_K, par.BehavWinner_TrajectoryMEAN);
par.BehavWinner_TrajectoryMEAN.InField              = 'exit_time';
[Jxy_WinnerMeanExT,Jxy_TabSoloMeanExT,Jxy_TabJointMeanExT]       = BehavWinner_Trajectory_MEAN(Jxy_S,Jxy_K, par.BehavWinner_TrajectoryMEAN);
par.BehavWinner_TrajectoryMEAN.InField              = 'ExitCorr';
[Jxy_WinnerMeanEC,Jxy_TabSoloMeanEC,Jxy_TabJointMeanEC]          = BehavWinner_Trajectory_MEAN(Jxy_S,Jxy_K, par.BehavWinner_TrajectoryMEAN);
par.BehavWinner_TrajectoryMEAN.InField              = 'CMT';
[Jxy_WinnerMeanCMT,Jxy_TabSoloMeanCMT,Jxy_TabJointMeanCMT]       = BehavWinner_Trajectory_MEAN(Jxy_S,Jxy_K, par.BehavWinner_TrajectoryMEAN);
% % Intra-monkey comparison tables / Within-monkey statistics
par.BehavWinner_Trajectory_MEAN_INTRA.num_dir                   = num_dir;
par.BehavWinner_Trajectory_MEAN_INTRA.InField                   = 'angular_error';
[Jxy_WinnerMeanAE_intra,Jxy_TabMeanAE_S,Jxy_TabMeanAE_K]        = BehavWinner_Trajectory_MEAN_INTRA(Jxy_S,Jxy_K, par.BehavWinner_Trajectory_MEAN_INTRA);
par.BehavWinner_Trajectory_MEAN_INTRA.InField                   = 'angular_error_magnitude';
[Jxy_WinnerMeanAE_intra_mod,Jxy_TabMeanAE_S_mod,Jxy_TabMeanAE_K_mod]  = BehavWinner_Trajectory_MEAN_INTRA(Jxy_S,Jxy_K, par.BehavWinner_Trajectory_MEAN_INTRA);
par.BehavWinner_Trajectory_MEAN_INTRA.InField                   = 'exit_time';
[Jxy_WinnerMeanExT_intra,Jxy_TabMeanExT_S,Jxy_TabMeanExT_K]     = BehavWinner_Trajectory_MEAN_INTRA(Jxy_S,Jxy_K, par.BehavWinner_Trajectory_MEAN_INTRA);
par.BehavWinner_Trajectory_MEAN_INTRA.InField                   = 'ExitCorr';
[Jxy_WinnerMeanEC_intra,Jxy_TabMeanEC_S,Jxy_TabMeanEC_K]        = BehavWinner_Trajectory_MEAN_INTRA(Jxy_S,Jxy_K, par.BehavWinner_Trajectory_MEAN_INTRA);
par.BehavWinner_Trajectory_MEAN_INTRA.InField                   = 'CMT';
[Jxy_WinnerMeanCMT_intra,Jxy_TabMeanCMT_S,Jxy_TabMeanCMT_K]     = BehavWinner_Trajectory_MEAN_INTRA(Jxy_S,Jxy_K, par.BehavWinner_Trajectory_MEAN_INTRA);

%% Other Params
par.BehavWinner.num_dir                         = num_dir;
par.BehavWinner.mode                            = 'time';
[RT_Winner,RT_TabSolo,RT_TabJoint]              = BehavWinner(RT_session, par.BehavWinner);

par.BehavWinner.mode                            = 'velocity';
[PV_Winner,PV_TabSolo,PV_TabJoint]              = BehavWinner(PV_session, par.BehavWinner);

par.BehavWinner.mode                            = 'time';
[PVT_Winner,PVT_TabSolo,PVT_TabJoint]           = BehavWinner(PVT_session, par.BehavWinner);

par.BehavWinner.mode                            = 'time';
[ST_Winner,ST_TabSolo,ST_TabJoint]              = BehavWinner(ST_session, par.BehavWinner);

%% MEAN FOR SESSION AND DIRECTION
par.meanSessionDir.SessName                 = session_list;
par.meanSessionDir.num_cond                 = num_cond;
par.meanSessionDir.num_dir                  = num_dir;
par.meanSessionDir.millisec                 = true;
[RT_SessionMEAN,RT_SessionSTD]              = meanSessionDir(RT_session,par.meanSessionDir);

par.meanSessionDir.SessName                 = session_list;
par.meanSessionDir.num_cond                 = num_cond;
par.meanSessionDir.num_dir                  = num_dir;
par.meanSessionDir.millisec                 = true;
[ST_sessionMEAN,ST_sessionSTD]              = meanSessionDir(ST_session,par.meanSessionDir);

par.meanSessionDir.SessName                 = session_list;
par.meanSessionDir.num_cond                 = num_cond;
par.meanSessionDir.num_dir                  = num_dir;
par.meanSessionDir.millisec                 = false;
[PV_sessionMEAN,PV_sessionSTD]              = meanSessionDir(PV_session,par.meanSessionDir);

par.meanSessionDir.SessName                 = session_list;
par.meanSessionDir.num_cond                 = num_cond;
par.meanSessionDir.num_dir                  = num_dir;
par.meanSessionDir.millisec                 = true;
[PVT_sessionMEAN,PVT_sessionSTD]            = meanSessionDir(PVT_session,par.meanSessionDir);
%% Mean Winner and Statistic
par.BehavWinnerMEAN.num_dir                         = num_dir;
par.BehavWinnerMEAN.mode                            = 'time';
[RT_WinnerMEAN,RT_TabSoloMEAN,RT_TabJointMEAN]      = BehavWinnerMEAN(RT_SessionMEAN, par.BehavWinnerMEAN);

par.BehavWinnerMEAN.mode                            = 'velocity';
[PV_WinnerMEAN,PV_TabSoloMEAN,PV_TabJointMEAN]      =  BehavWinnerMEAN(PV_sessionMEAN, par.BehavWinnerMEAN);

par.BehavWinnerMEAN.mode                            = 'time';
[PVT_WinnerMEAN,PVT_TabSoloMEAN,PVT_TabJointMEAN]   = BehavWinnerMEAN(PVT_sessionMEAN, par.BehavWinnerMEAN);

par.BehavWinnerMEAN.mode                            = 'time';
[ST_WinnerMEAN,ST_TabSoloMEAN,ST_TabJointMEAN]      = BehavWinnerMEAN(ST_sessionMEAN, par.BehavWinnerMEAN);
% intra
par.BehavWinnerMEAN_INTRA.num_dir                   = num_dir;
par.BehavWinnerMEAN_INTRA.mode                      = 'time';
[RT_WinnerMEAN_Intra,RT_TabMEAN_S,RT_TabMEAN_K]           = BehavWinnerMEAN_INTRA(RT_SessionMEAN, par.BehavWinnerMEAN_INTRA);

par.BehavWinnerMEAN_INTRA.mode                      = 'velocity';
[PV_WinnerMEAN_Intra,PV_TabMEAN_S,PV_TabMEAN_K]           = BehavWinnerMEAN_INTRA(PV_sessionMEAN, par.BehavWinnerMEAN_INTRA);

par.BehavWinnerMEAN_INTRA.mode                      = 'time';
[PVT_WinnerMEAN_Intra,PVT_TabMEAN_S,PVT_TabMEAN_K]        = BehavWinnerMEAN_INTRA(PVT_sessionMEAN, par.BehavWinnerMEAN_INTRA);

par.BehavWinnerMEAN_INTRA.mode                      = 'time';
[ST_WinnerMEAN_Intra,ST_TabMEAN_S,ST_TabMEAN_K]        = BehavWinnerMEAN_INTRA(ST_sessionMEAN, par.BehavWinnerMEAN_INTRA);

%% Arrange matrix
par.matrixArrange.num_cond  = num_cond;
par.matrixArrange.num_dir   = num_dir;
RT_total                    = matrixArrange(RT_SessionMEAN,par.matrixArrange);
ST_total                    = matrixArrange(ST_sessionMEAN,par.matrixArrange);
PV_total                    = matrixArrange(PV_sessionMEAN,par.matrixArrange);
PVT_total                   = matrixArrange(PVT_sessionMEAN,par.matrixArrange);

% total mean
par.meanOutEval.num_cond                                 = num_cond;
par.meanOutEval.num_dir                                  = num_dir;
[RT_total_Mean,~,RT_total_SE,RTplot_maxY]                = meanOutEval(RT_total,par.meanOutEval);
[ST_total_Mean,~,ST_total_SE,STplot_maxY]                = meanOutEval(ST_total,par.meanOutEval);
[PV_total_Mean,~,PV_total_SE,PV_total_maxY]              = meanOutEval(PV_total,par.meanOutEval);
[PVT_total_Mean,~,PVT_total_SE,PVT_total_maxY]           = meanOutEval(PVT_total,par.meanOutEval);

%% Build output structure
Data.Behav.Jxy_S        = Jxy_S;
Data.Behav.Jxy_K        = Jxy_K;
Data.Behav.RT_total     = RT_total;
Data.Behav.PV_total     = PV_total;
Data.Behav.PVT_total    = PVT_total;
Data.Behav.ST_total     = ST_total;
Data.Behav.Jxy_session  = Jxy_session;
Data.Behav.RT_session   = RT_session;
Data.Behav.PV_session   = PV_session;
Data.Behav.PVT_session  = PVT_session;
Data.Behav.ST_session   = ST_session;
%% Tab
Data.Behav.Jxy_TabSoloAE    = Jxy_TabSoloAE;
Data.Behav.Jxy_TabJointAE   = Jxy_TabJointAE;
Data.Behav.Jxy_TabSoloAE_abs    = Jxy_TabSoloAE_abs;
Data.Behav.Jxy_TabJointAE_abs   = Jxy_TabJointAE_abs;
Data.Behav.Jxy_TabSoloExT   = Jxy_TabSoloExT;
Data.Behav.Jxy_TabJointExT  = Jxy_TabJointExT;
Data.Behav.Jxy_TabSoloEC    = Jxy_TabSoloEC;
Data.Behav.Jxy_TabJointEC   = Jxy_TabJointEC;
Data.Behav.Jxy_TabSoloCMT   = Jxy_TabSoloCMT;
Data.Behav.Jxy_TabJointCMT  = Jxy_TabJointCMT;
Data.Behav.RT_TabSolo       = RT_TabSolo;
Data.Behav.RT_TabJoint      = RT_TabJoint;
Data.Behav.PV_TabSolo       = PV_TabSolo;
Data.Behav.PV_TabJoint      = PV_TabJoint;
Data.Behav.PVT_TabSolo      = PVT_TabSolo;
Data.Behav.PVT_TabJoint     = PVT_TabJoint;
Data.Behav.ST_TabSolo       = ST_TabSolo;
Data.Behav.ST_TabJoint      = ST_TabJoint;
 
% mean
Data.Behav.Jxy_TabSoloMeanAE    = Jxy_TabSoloMeanAE;
Data.Behav.Jxy_TabJointMeanAE   = Jxy_TabJointMeanAE;
Data.Behav.Jxy_TabSoloMeanAE_abs    = Jxy_TabSoloMeanAE_abs;
Data.Behav.Jxy_TabJointMeanAE_abs   = Jxy_TabJointMeanAE_abs;
Data.Behav.Jxy_TabSoloMeanExT   = Jxy_TabSoloMeanExT;
Data.Behav.Jxy_TabJointMeanExT  = Jxy_TabJointMeanExT;
Data.Behav.Jxy_TabSoloMeanEC    = Jxy_TabSoloMeanEC;
Data.Behav.Jxy_TabJointMeanEC   = Jxy_TabJointMeanEC;
Data.Behav.Jxy_TabSoloMeanCMT   = Jxy_TabSoloMeanCMT;
Data.Behav.Jxy_TabJointMeanCMT  = Jxy_TabJointMeanCMT;
Data.Behav.RT_TabSoloMEAN       = RT_TabSoloMEAN;
Data.Behav.RT_TabJointMEAN      = RT_TabJointMEAN;
Data.Behav.PV_TabSoloMEAN       = PV_TabSoloMEAN;
Data.Behav.PV_TabJointMEAN      = PV_TabJointMEAN;
Data.Behav.PVT_TabSoloMEAN      = PVT_TabSoloMEAN;
Data.Behav.PVT_TabJointMEAN     = PVT_TabJointMEAN;
Data.Behav.ST_TabSoloMEAN       = ST_TabSoloMEAN;
Data.Behav.ST_TabJointMEAN      = ST_TabJointMEAN;

Data.Behav.Jxy_TabMeanAE_S      = Jxy_TabMeanAE_S;
Data.Behav.Jxy_TabMeanAE_K      = Jxy_TabMeanAE_K;
Data.Behav.Jxy_TabMeanAE_S_mod      = Jxy_TabMeanAE_S_mod;
Data.Behav.Jxy_TabMeanAE_K_mod      = Jxy_TabMeanAE_K_mod;
Data.Behav.Jxy_TabMeanExT_S     = Jxy_TabMeanExT_S;
Data.Behav.Jxy_TabMeanExT_K     = Jxy_TabMeanExT_K;
Data.Behav.Jxy_TabMeanEC_S      = Jxy_TabMeanEC_S;
Data.Behav.Jxy_TabMeanEC_K      = Jxy_TabMeanEC_K;
Data.Behav.Jxy_TabMeanCMT_S     = Jxy_TabMeanCMT_S;
Data.Behav.Jxy_TabMeanCMT_K     = Jxy_TabMeanCMT_K;
Data.Behav.RT_TabMEAN_S         = RT_TabMEAN_S;
Data.Behav.RT_TabMEAN_K         = RT_TabMEAN_K;
Data.Behav.PV_TabMEAN_S         = PV_TabMEAN_S;
Data.Behav.PV_TabMEAN_K         = PV_TabMEAN_K;
Data.Behav.PVT_TabMEAN_S        = PVT_TabMEAN_S;
Data.Behav.PVT_TabMEAN_K        = PVT_TabMEAN_K;
Data.Behav.ST_TabMEAN_S         = ST_TabMEAN_S;
Data.Behav.ST_TabMEAN_K         = ST_TabMEAN_K;


%% SAVE DATA

output_dir = fullfile(pwd,directory_save);

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

save_dir = fullfile(output_dir,file_name);

%% SAVE DATA STRUCTURE

save(save_dir,'Data','-v7.3');

fprintf('\n====================================================\n');
fprintf('Data successfully saved:\n');
fprintf('%s\n',save_dir);
fprintf('====================================================\n');