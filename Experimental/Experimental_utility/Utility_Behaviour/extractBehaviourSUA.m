function behaviour = extractBehaviourSUA(session_name,pms)
% EXTRACTBEHAVIORSUA Extract behavioral metrics from a SUA session.
%
% This function loads a SUA recording session, extracts kinematic and
% behavioral variables for each trial, organizes trials by condition and
% movement direction, and computes summary behavioral metrics.
%
% INPUTS
% -------
% session_name : char | string
%       Path to the SUA session file (.mat).
%
% pms : struct (optional)
%       Parameter structure.
%
%       interval       : [start end] analysis window (s)
%                        Default: [0 0.5]
%
%       AlignEvent     : event code used for alignment
%                        Default: 2024
%
%       LockEvent      : event code used for temporal locking
%                        Default: 2024
%
%       tstartExtract  : extraction start time relative to AlignEvent (s)
%                        Default: -0.7
%
%       tendExtract    : extraction end time relative to AlignEvent (s)
%                        Default: 5.0
%
% OUTPUT
% -------
% behavior : struct
%
%       behavior.S
%           Metrics for monkey S
%
%       behavior.K
%           Metrics for monkey K
%
%       Each contains 24 entries:
%           3 conditions × 8 directions
%
%       Fields:
%           trialType
%           Direction
%           Condition
%           trialName
%           Jxy
%           RT
%           PV
%           PVT
%           AE
%           ExT
%           EC
%           CMT
%           ICDmax
%           ICDmean
%           ICDauc
% Author: Mirco Frosolone
% Created: 2026-07-16
% DEPENDENCIES
% -------------------------------------------------------------------------
% getClassInfo
% RTextract
% ETextract
% alignBattaglia
% findSKdirection
% TrajXYextrat4eval
% XYtrajectoryEval
% RTextrat4eval
% PVExtrat4eval
% PVTExtrat4eval
% extractVarParams
% computeTrajectoryExitDir
% JxyMatrix
% extractBehavDeltas
% extractBehavDeltasPlus
% build24BehaviourStruct

defaults.interval      = [0 0.5];
defaults.AlignEvent    = 2024;
defaults.LockEvent     = 2024;
defaults.tstartExtract = -0.7;
defaults.tendExtract   = 5.0;

if nargin < 2
    pms = struct();
end

fields = fieldnames(defaults);

for k = 1:numel(fields)
    if ~isfield(pms,fields{k}) || isempty(pms.(fields{k}))
        pms.(fields{k}) = defaults.(fields{k});
    end
end

interval      = pms.interval;
AlignEvent    = pms.AlignEvent;
LockEvent     = pms.LockEvent;
tstartExtract = pms.tstartExtract;
tendExtract   = pms.tendExtract;

validateattributes(session_name,...
    {'char','string'},...
    {'nonempty'},...
    mfilename,...
    'session_name');

if ~exist([session_name,'.mat'],'file')
    error('Session file not found:\n%s',session_name);
end

try
    d = load([session_name,'.mat']);
catch ME
    error('Unable to load session file "%s":\n%s',...
        session_name,...
        ME.message);
end

if ~isfield(d,'Trials')
    error('File "%s" does not contain a ''Trials'' variable.',session_name);
end

Trials  = d.Trials;

[~, ~, klabels, success] = getClassInfo(Trials);
data_SUA     = Trials(success);
klabels      = klabels(success,:);

idMatch = regexp(session_name,'\d+','match','once');

if isempty(idMatch)
    warning('No session ID detected in file name.');
    sessionID = NaN;
else
    sessionID = str2double(idMatch);
end

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

% Positional offsets used to recenter trajectories:
% [Jx Jy Ex Ey Kx Ky Fx Fy]
shiftRawSUA = zeros(8,1);
% Monkey S offset
shiftRawSUA(1:2) = [0.5; -3.8];
% Monkey K offset
shiftRawSUA(5:6) = [1.0; -3.0];

for iTrial=1:numel(data_SUA)
    data_SUA(iTrial).trialId = iTrial;
    data_SUA(iTrial).trialName        = Labels{data_SUA(iTrial).trialType};
    data_SUA(iTrial).JXYEXY = data_SUA(iTrial).JXYEXY - shiftRawSUA;
end
%% Step 1: Arrange data
data_trials = findSKdirection(data_SUA);

num_dir         = 8; % number of direction
num_cond        = 3; % number of conditions: 1-SoloS 2-SoloK 3-Joint S-K

% ICD
[ICD_trials, ~] = ICD_trace(data_trials);
ICD_all         = ICD_statsEval(ICD_trials);

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
        fieldName = sprintf('dir%d',ndir);
        Move(cd).(fieldName) = data_trials(Condition{cd,ndir});
        lenM(cd,ndir)                           = length(data_trials(Condition{cd,ndir}));
    end
end

%% Trajectory Error
par.TrajXYextrat4eval.num_cond          = num_cond;
par.TrajXYextrat4eval.num_dir           = num_dir;
par.TrajXYextrat4eval.fieldname         = 'JXYEXY';
par.TrajXYextrat4eval.time              = 'timeJXYEXY';
[JxyS_trials,JxyK_trials,Jxy_trialsname,Jxy_trialstime,RT_Jxy_trials]    = TrajXYextrat4eval(Move,par.TrajXYextrat4eval);
[Jxy_trajS_trials,Jxy_trajK_trials]                   = XYtrajectoryEval(JxyS_trials,JxyK_trials);
Jxy_trials.session             = session_name;
Jxy_trials.S                   = JxyS_trials;
Jxy_trials.K                   = JxyK_trials;
Jxy_trials.name                = Jxy_trialsname;
Jxy_trials.time                = Jxy_trialstime;
Jxy_trials.RT_Jxy              = RT_Jxy_trials;
Jxy_trials.TrajS               = Jxy_trajS_trials;
Jxy_trials.TrajK               = Jxy_trajK_trials;

par.RTextrat4eval.num_cond  = num_cond;
par.RTextrat4eval.num_dir   = num_dir;
par.RTextrat4eval.idsess    = sessionID;
par.RTextrat4eval.fieldname = 'RT';
[RT_trials,~,~,~,~]     = RTextrat4eval(Move,lenM,par.RTextrat4eval);
par.extractVarParams.time = 'ms';
RT = extractVarParams(RT_trials,par.extractVarParams);

par.PVExtrat4eval.num_cond       = num_cond;
par.PVExtrat4eval.num_dir        = num_dir;
par.PVExtrat4eval.idsess         = sessionID;
par.PVExtrat4eval.fieldname      = 'JXYEXY';
par.PVExtrat4eval.time           = 'timeJXYEXY';
par.PVExtrat4eval.windowSize     = 50;
[PV_trials,~,~,~,~]                = PVExtrat4eval(Move,lenM,par.PVExtrat4eval);
par.extractVarParams.time = 'NaN';
PV = extractVarParams(PV_trials,par.extractVarParams);

par.PVTExtrat4eval.num_cond            = num_cond;
par.PVTExtrat4eval.num_dir             = num_dir;
par.PVTExtrat4eval.idsess              = sessionID;
par.PVTExtrat4eval.fieldname           = 'JXYEXY';
par.PVTExtrat4eval.time                = 'timeJXYEXY';
par.PVTExtrat4eval.windowSize          = 50;
[PVT_trials,~,~,~,~]                   = PVTExtrat4eval(Move,lenM,par.PVTExtrat4eval);
par.extractVarParams.time = 'ms';
PVT = extractVarParams(PVT_trials,par.extractVarParams);

par.computeTrajectoryExitDir.r_dim = 1.81;
Jxy_Result = computeTrajectoryExitDir(Jxy_trials,par.computeTrajectoryExitDir);
[Jxy_S, Jxy_K] = JxyMatrix(Jxy_Result);
[~, ~, AE] = extractBehavDeltas(Jxy_S, Jxy_K, 'angular_error');
[~, ~, AE_abs] = extractBehavDeltas(Jxy_S, Jxy_K, 'angular_error_magnitude');
[~, ~, CMT] = extractBehavDeltas(Jxy_S, Jxy_K, 'CMT');
[~, ~, EC] = extractBehavDeltasPlus(Jxy_S, Jxy_K, 'ExitCorr');
[~, ~, ExT] = extractBehavDeltas(Jxy_S, Jxy_K, 'exit_time');

%% Build output structure
% Each call appends a specific behavioral metric to the
% 24 condition-direction entries.
behaviour = build24BehaviourStruct(AE);
behaviour = build24BehaviourStruct(AE_abs,behaviour);
behaviour = build24BehaviourStruct(RT, behaviour);
behaviour = build24BehaviourStruct(PV, behaviour);
behaviour = build24BehaviourStruct(PVT, behaviour);
behaviour = build24BehaviourStruct(CMT, behaviour);
behaviour = build24BehaviourStruct(EC, behaviour);
behaviour = build24BehaviourStruct(ExT, behaviour);
% ICD
behaviour.ICD = build24BehaviourICDStruct(ICD_all);

