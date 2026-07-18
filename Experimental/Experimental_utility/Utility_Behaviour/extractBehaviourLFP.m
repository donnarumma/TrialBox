function behaviour = extractBehaviourLFP(session_name,pms)
% EXTRACTBEHAVIOURLFP Extract behavioral metrics from a LFP session.
%
% This function loads an LFP recording session, extracts kinematic and
% behavioral variables for each trial, organizes trials by condition and
% movement direction, and computes summary behavioral metrics.
%
% INPUTS
% -------
% session_name : char | string
%       Path to the LFP session file.
%
% pms : struct (optional)
%
%       interval       : [start end] analysis window (s)
%                        Default: [0 0.5]
%
%       AlignEvent     : event code used for alignment
%                        Default: 2024
%
%       LockEvent      : event code used for locking
%                        Default: 2024
%
%       selectionLFPs  : channel selection method
%                        Default: 'maximum'
%
% OUTPUT
% -------
% behaviour : struct
%
%       behaviour.S
%       behaviour.K
%
%       Each contains 24 entries:
%           3 conditions × 8 directions
%
% Author: Mirco Frosolone
% Created: 2026-07-10
% DEPENDENCIES
% -------------------------------------------------------------------------
% selectionLFPs
% BattagliaArrangeTrialsParams
% BattagliaArrangeAlignTrials
% alignBattaglia
% trimAlignedField
% resampleJXYtoLFP
% getJointMonkeysLabels
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
defaults.selectionLFPs = 'maximum';

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

tStart                  = interval(1);
tEnd                    = interval(2);

if isfield(pms, 'selectionLFPs') && ~isempty(pms.selectionLFPs)
        lfpSelection.method = pms.selectionLFPs;
    else
        lfpSelection.method = 'maximum';
end

validateattributes(session_name,...
    {'char','string'},...
    {'nonempty'},...
    mfilename,...
    'session_name');

[selS, selK, ~, ~] = selectionLFPs(session_name, lfpSelection);

if isempty(selS) && isempty(selK)
    error('Nessun LFP selezionato per la sessione: %s', session_name);
end

idMatch = regexp(session_name,'\d+','match','once');

if isempty(idMatch)
    warning('No session ID detected in file name.');
    sessionID = NaN;
else
    sessionID = str2double(idMatch);
end
%% Step 0: Arrange Trials
fprintf('Step 0: arrange trials\n');
par.BattagliaArrangeTrials              = BattagliaArrangeTrialsParams();
par.BattagliaArrangeTrials.selS         = selS;
par.BattagliaArrangeTrials.selK         = selK;
par.BattagliaArrangeTrials.session_name = session_name;
par.BattagliaArrangeTrials.tstart       = -0.7;
par.BattagliaArrangeTrials.tstop        = 2.2;
par.BattagliaArrangeTrials.AlignEvent   = pms.AlignEvent;
par.BattagliaArrangeTrials.InField      = 'LFP';

data_trials = BattagliaArrangeAlignTrials(par.BattagliaArrangeTrials);

% Delete label 0 trials
idx_empty = find(arrayfun(@(x) isempty(x.trialType), data_trials));

if ~isempty(idx_empty)
    par.badTrialsInfo.name     = session_name;
    par.badTrialsInfo.chamber  = {data_trials(idx_empty).Chamber};
    par.badTrialsInfo.trialsId = [data_trials(idx_empty).trialId];
    par.badTrialsInfo.trials   = data_trials(idx_empty);
else
    par.badTrialsInfo.name     = session_name;
    par.badTrialsInfo.chamber  = data_trials(1).Chamber;
    par.badTrialsInfo.trialsId = [];
    par.badTrialsInfo.trials   = [];
end

data_trials(idx_empty) = [];

%% Align LFP signal
par.alignBattaglia.AlignEvent = AlignEvent;
par.alignBattaglia.LockEvent  = LockEvent;
par.alignBattaglia.InField    = 'LFP'; % e.g. 'LFP', 'JXYEXY'
par.alignBattaglia.OutField   = 'LFP'; % e.g. 'LFP', 'JXYEXY'
data_trials = alignBattaglia(data_trials, par.alignBattaglia);

%% Time trimming after alignment
par.trimAlignedField.InField = 'LFP'; % e.g. 'LFP', 'JXYEXY'
par.trimAlignedField.tStart  = tStart;
par.trimAlignedField.tEnd    = tEnd;
data_trials = trimAlignedField(data_trials, par.trimAlignedField);

%% Align JXYEXY signal
par.alignBattaglia.AlignEvent = AlignEvent;
par.alignBattaglia.LockEvent  = LockEvent;
par.alignBattaglia.InField    = 'JXYEXY'; % e.g. 'LFP', 'JXYEXY'
par.alignBattaglia.OutField   = 'JXYEXY'; % e.g. 'LFP', 'JXYEXY'
data_trials = alignBattaglia(data_trials, par.alignBattaglia);

%% Time trimming after alignment
par.trimAlignedField.InField = 'JXYEXY'; % e.g. 'LFP', 'JXYEXY'
par.trimAlignedField.tStart  = tStart;
par.trimAlignedField.tEnd    = tEnd;
data_trials = trimAlignedField(data_trials, par.trimAlignedField);
% add trialName
[~,Labels] = getJointMonkeysLabels(1:24);
for iTrial=1:numel(data_trials)
    data_trials(iTrial).trialName        = Labels{data_trials(iTrial).trialType};
end
%% resample JXY
data_trials = resampleJXYtoLFP(data_trials);
data_trials = getSKconditiondirection(data_trials);
%% Step 1: Arrange data
num_dir         = 8; % number of direction
num_cond        = 3; % number of conditions: 1-SoloS 2-SoloK 3-Joint S-K

data_trials = findSKdirection(data_trials);

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
