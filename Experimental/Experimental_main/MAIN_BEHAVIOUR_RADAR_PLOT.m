%% MAIN_BEHAVIOUR_RADAR_PLOTS
% MAIN_BEHAVIOUR_RADAR_PLOTS.m
%
% Generates radar plots for behavioural metrics obtained from
% Behaviour_*_TotalSessionsData MAT files.
%
% Generated figures:
%   - Solo S vs Solo K
%   - Joint S vs Joint K
%   - Within-monkey comparisons
%   - Four-condition comparisons
%   - ICD summaries
%
% Supported chambers:
%   - Frontal
%   - Parietal
%
% Author: Mirco Frosolone
% Date: 2026-07-06
%
% Metric abbreviations:
% AE   = Angular Error
% ExT  = Exit Time
% EC   = Exit Correction
% CMT  = Corrective Movement Time
% RT   = Reaction Time
% PV   = Peak Velocity
% PVT  = Peak Velocity Time
% ICD  = Inter-Cursor Distance

clear;
close all;

%% Behavioural Radar Plots
origState = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','off');

params = struct();
params.chamber = 'Frontal';

behavFile = sprintf( ...
    'Behaviour_%s_TotalSessionsData_0_05s.mat', ...
    upper(params.chamber));

icdFile = sprintf( ...
    'ICD_%s_MovOn_200.mat', ...
    upper(params.chamber));
switch lower(params.chamber)

    case 'frontal'
        areaFolder = 'FRONTAL';

    case 'parietal'
        areaFolder = 'PARIETAL';

    otherwise
        error('Unknown chamber: %s', params.chamber);
end

save_dir = fullfile( ...
    pwd,...
    'BEHAVIOR',...
    't0_t500',...
    areaFolder,...
    'POLARPLOT');

S   = load(behavFile);
Data = S.Data;

icd = load(icdFile);

%% Palette
S = load('Palette.mat');
Palette = S.Palette;

FS_TITLE  = Palette.fontStyle.font_title;
FS_LEGEND = Palette.fontStyle.font_legend;
FS_TICKS  = Palette.fontStyle.font_ticks;

LS_ACT    = Palette.lineStyle.Act;
LS_JOINT  = Palette.lineStyle.Joint;

%% Data extraction
Jxy_S       = Data.Behav.Jxy_S;
Jxy_K       = Data.Behav.Jxy_K;
RT_total    = Data.Behav.RT_total;
PV_total    = Data.Behav.PV_total;
PVT_total   = Data.Behav.PVT_total;

frac  = 0.15;  % 15% of max

%% AE
[~, ~, AE] = extractBehavDeltas(Jxy_S, Jxy_K, 'angular_error');
AE_SoloK  = mean(AE.K.Solo,  'omitnan')';
AE_SoloS  = mean(AE.S.Solo,  'omitnan')';
AE_JointK = mean(AE.K.Joint, 'omitnan')';
AE_JointS = mean(AE.S.Joint, 'omitnan')';

maxAE = max([AE_SoloS;
             AE_SoloK;
             AE_JointS;
             AE_JointK], [], 'omitnan');
deltaAE     = frac * maxAE;
maxAE       = maxAE + deltaAE;
minAE = min([AE_SoloS;
             AE_SoloK;
             AE_JointS;
             AE_JointK], [], 'omitnan');
minAE       = minAE - deltaAE;
%% AE_abs
[~, ~, AE_abs] = extractBehavDeltas(Jxy_S, Jxy_K, 'angular_error_magnitude');
AE_abs_SoloK  = mean(abs(AE_abs.K.Solo),  'omitnan')';
AE_abs_SoloS  = mean(abs(AE_abs.S.Solo),  'omitnan')';
AE_abs_JointK = mean(abs(AE_abs.K.Joint), 'omitnan')';
AE_abs_JointS = mean(abs(AE_abs.S.Joint), 'omitnan')';

maxAE_abs = max([AE_abs_SoloS;
             AE_abs_SoloK;
             AE_abs_JointS;
             AE_abs_JointK], [], 'omitnan');
deltaAE_abs     = frac * maxAE_abs;
maxAE_abs       = maxAE_abs + deltaAE_abs;
minAE_abs = min([AE_abs_SoloS;
             AE_abs_SoloK;
             AE_abs_JointS;
             AE_abs_JointK], [], 'omitnan');
minAE_abs       = minAE_abs - deltaAE_abs;
%% ExT
[~, ~, ExT] = extractBehavDeltas(Jxy_S, Jxy_K, 'exit_time');
ExT_SoloK   = 1000*mean(ExT.K.Solo, 'omitnan')';
ExT_SoloS   = 1000*mean(ExT.S.Solo, 'omitnan')';
ExT_JointK  = 1000*mean(ExT.K.Joint, 'omitnan')';
ExT_JointS  = 1000*mean(ExT.S.Joint, 'omitnan')';

maxExT      = max([ExT_SoloS; ExT_SoloK; ExT_JointS; ExT_JointK], [], 'omitnan');
deltaExT    = frac * maxExT;
maxExT      = maxExT + deltaExT;
minExT      = min([ExT_SoloS; ExT_SoloK; ExT_JointS; ExT_JointK], [], 'omitnan');
minExT      = minExT - deltaExT;

%% EC
[~, ~, EC] = extractBehavDeltasPlus(Jxy_S, Jxy_K, 'ExitCorr');
EC_SoloK    = mean(EC.K.Solo, 'omitnan')';
EC_SoloS    = mean(EC.S.Solo, 'omitnan')';
EC_JointK   = mean(EC.K.Joint, 'omitnan')';
EC_JointS   = mean(EC.S.Joint, 'omitnan')';

maxEC       = max([EC_SoloS; EC_SoloK; EC_JointS; EC_JointK], [], 'omitnan');
deltaEC     = frac * maxEC;
maxEC       = maxEC + deltaEC;
minEC       = min([EC_SoloS; EC_SoloK; EC_JointS; EC_JointK], [], 'omitnan');
minEC       = minEC - deltaEC;

%% CMT
[~, ~, CMT] = extractBehavDeltas(Jxy_S, Jxy_K, 'CMT');
CMT_SoloK   = 1000*mean(CMT.K.Solo, 'omitnan')';
CMT_SoloS   = 1000*mean(CMT.S.Solo, 'omitnan')';
CMT_JointK  = 1000*mean(CMT.K.Joint, 'omitnan')';
CMT_JointS  = 1000*mean(CMT.S.Joint, 'omitnan')';

maxCMT      = max([CMT_SoloS; CMT_SoloK; CMT_JointS; CMT_JointK], [], 'omitnan');
deltaCMT    = frac * maxCMT;
maxCMT      = maxCMT + deltaCMT;
minCMT      = min([CMT_SoloS; CMT_SoloK; CMT_JointS; CMT_JointK], [], 'omitnan');
minCMT      = minCMT - deltaCMT;

%% RT
[~, ~, RT] = extractFromTotalStruct(RT_total, 'less');
RT_SoloK    = mean(RT.K.Solo, 'omitnan')';
RT_SoloS    = mean(RT.S.Solo, 'omitnan')';
RT_JointK   = mean(RT.K.Joint, 'omitnan')';
RT_JointS   = mean(RT.S.Joint, 'omitnan')';

maxRT       = max([RT_SoloS; RT_SoloK; RT_JointS; RT_JointK], [], 'omitnan');
deltaRT     = frac * maxRT;
maxRT       = maxRT + deltaRT;
minRT       = min([RT_SoloS; RT_SoloK; RT_JointS; RT_JointK], [], 'omitnan');
minRT       = minRT - deltaRT;

%% PV
[~, ~, PV] = extractFromTotalStruct(PV_total, 'more');
PV_SoloK    = mean(PV.K.Solo, 'omitnan')';
PV_SoloS    = mean(PV.S.Solo, 'omitnan')';
PV_JointK   = mean(PV.K.Joint, 'omitnan')';
PV_JointS   = mean(PV.S.Joint, 'omitnan')';

maxPV       = max([PV_SoloS; PV_SoloK; PV_JointS; PV_JointK], [], 'omitnan');
deltaPV     = frac * maxPV;
maxPV       = maxPV + deltaPV;
minPV       = min([PV_SoloS; PV_SoloK; PV_JointS; PV_JointK], [], 'omitnan');
minPV       = minPV - deltaPV;

%% PVT
[~, ~, PVT] = extractFromTotalStruct(PVT_total, 'less');
PVT_SoloK   = mean(PVT.K.Solo, 'omitnan')';
PVT_SoloS   = mean(PVT.S.Solo, 'omitnan')';
PVT_JointK  = mean(PVT.K.Joint, 'omitnan')';
PVT_JointS  = mean(PVT.S.Joint, 'omitnan')';

maxPVT      = max([PVT_SoloS; PVT_SoloK; PVT_JointS; PVT_JointK], [], 'omitnan');
deltaPVT    = frac * maxPVT;
maxPVT      = maxPVT + deltaPVT;
minPVT      = min([PVT_SoloS; PVT_SoloK; PVT_JointS; PVT_JointK], [], 'omitnan');
minPVT      = minPVT - deltaPVT;

%% Common save dir
rootSaveDir = save_dir;

ttestSaveDir = fullfile(rootSaveDir,'T-Test');

SOLO_FIG = 'ActS_vs_ActK';
JOINT_FIG = 'JointS_vs_JointK';
commonRadPlot = struct();

commonRadPlot.font_title  = FS_TITLE;
commonRadPlot.font_legend = FS_LEGEND;
commonRadPlot.font_ticks  = FS_TICKS;
commonRadPlot.mod         = 'NaN';

par.RadPlot = commonRadPlot;
par.RadPlotOwn = commonRadPlot;
par.RadMultiplePlot = commonRadPlot;
%% =========================================================================
%% T-TEST: SOLO S vs SOLO K
%% =========================================================================
% ========================= AE =========================
par.RadPlot.save_dir = fullfile(ttestSaveDir,'AE','Single');
par.RadPlot.unitMeasure  = '°';
par.RadPlot.figureName   = SOLO_FIG;
par.RadPlot.cond_names1  = {'Act S'};
par.RadPlot.cond_names2  = {'Act K'};
par.RadPlot.InField      = 'AE';
par.RadPlot.color1       = Palette.color.ActS;
par.RadPlot.color2       = Palette.color.ActK;
par.RadPlot.lineStyle1   = LS_ACT;
par.RadPlot.lineStyle2   = LS_ACT;
par.RadPlot.legendNames  = {{'Monkey S','AE','Act S'}, ...
                            {'Monkey K','AE','Act K'}};
par.RadPlot.p_vals       = Data.Behav.Jxy_TabSoloMeanAE.pvalue';
par.RadPlot.maxValue     = maxAE;
par.RadPlot.minValue     = minAE;
RadPlot(AE_SoloS', AE_SoloK', par.RadPlot)
% ========================= AE_abs =========================
par.RadPlot.save_dir = fullfile(ttestSaveDir,'AE_abs','Single');
par.RadPlot.unitMeasure  = '°';
par.RadPlot.figureName   = SOLO_FIG;
par.RadPlot.cond_names1  = {'Act S'};
par.RadPlot.cond_names2  = {'Act K'};
par.RadPlot.InField      = 'AE';
par.RadPlot.color1       = Palette.color.ActS;
par.RadPlot.color2       = Palette.color.ActK;
par.RadPlot.lineStyle1   = LS_ACT;
par.RadPlot.lineStyle2   = LS_ACT;
par.RadPlot.legendNames  = {{'Monkey S','AE_abs','Act S'}, ...
                            {'Monkey K','AE_abs','Act K'}};
par.RadPlot.p_vals       = Data.Behav.Jxy_TabSoloMeanAE_abs.pvalue';
par.RadPlot.maxValue     = maxAE_abs;
par.RadPlot.minValue     = minAE_abs;
RadPlot(AE_abs_SoloS', AE_abs_SoloK', par.RadPlot)
% ========================= ExT =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'ExT','Single');
par.RadPlot.unitMeasure  = 'ms';
par.RadPlot.figureName   = SOLO_FIG;
par.RadPlot.cond_names1  = {'Act S'};
par.RadPlot.cond_names2  = {'Act K'};
par.RadPlot.InField      = 'ExT';
par.RadPlot.color1       = Palette.color.ActS;
par.RadPlot.color2       = Palette.color.ActK;
par.RadPlot.lineStyle1   = LS_ACT;
par.RadPlot.lineStyle2   = LS_ACT;
par.RadPlot.legendNames  = {{'Monkey S','ExT','Act S'}, ...
                            {'Monkey K','ExT','Act K'}};
par.RadPlot.p_vals       = Data.Behav.Jxy_TabSoloMeanExT.pvalue';
par.RadPlot.maxValue     = maxExT;
par.RadPlot.minValue     = minExT;
RadPlot(ExT_SoloS', ExT_SoloK', par.RadPlot)

% ========================= EC =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'EC','Single');
par.RadPlot.unitMeasure  = '%';
par.RadPlot.figureName   = SOLO_FIG;
par.RadPlot.cond_names1  = {'Act S'};
par.RadPlot.cond_names2  = {'Act K'};
par.RadPlot.InField      = 'EC';
par.RadPlot.color1       = Palette.color.ActS;
par.RadPlot.color2       = Palette.color.ActK;
par.RadPlot.lineStyle1   = LS_ACT;
par.RadPlot.lineStyle2   = LS_ACT;
par.RadPlot.legendNames  = {{'Monkey S','EC','Act S'}, ...
                            {'Monkey K','EC','Act K'}};
par.RadPlot.p_vals       = Data.Behav.Jxy_TabSoloMeanEC.pvalue';
par.RadPlot.maxValue     = maxEC;
par.RadPlot.minValue     = minEC;
RadPlot(EC_SoloS', EC_SoloK', par.RadPlot)

% ========================= CMT =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'CMT','Single');
par.RadPlot.unitMeasure  = 'ms';
par.RadPlot.figureName   = SOLO_FIG;
par.RadPlot.cond_names1  = {'Act S'};
par.RadPlot.cond_names2  = {'Act K'};
par.RadPlot.InField      = 'CMT';
par.RadPlot.color1       = Palette.color.ActS;
par.RadPlot.color2       = Palette.color.ActK;
par.RadPlot.lineStyle1   = LS_ACT;
par.RadPlot.lineStyle2   = LS_ACT;
par.RadPlot.legendNames  = {{'Monkey S','CMT','Act S'}, ...
                            {'Monkey K','CMT','Act K'}};
par.RadPlot.p_vals       = Data.Behav.Jxy_TabSoloMeanCMT.pvalue';
par.RadPlot.maxValue     = maxCMT;
par.RadPlot.minValue     = minCMT;
RadPlot(CMT_SoloS', CMT_SoloK', par.RadPlot)

% ========================= RT =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'RT','Single');
par.RadPlot.unitMeasure  = 'ms';
par.RadPlot.figureName   = SOLO_FIG;
par.RadPlot.cond_names1  = {'Act S'};
par.RadPlot.cond_names2  = {'Act K'};
par.RadPlot.InField      = 'RT';
par.RadPlot.color1       = Palette.color.ActS;
par.RadPlot.color2       = Palette.color.ActK;
par.RadPlot.lineStyle1   = LS_ACT;
par.RadPlot.lineStyle2   = LS_ACT;
par.RadPlot.legendNames  = {{'Monkey S','RT','Act S'}, ...
                            {'Monkey K','RT','Act K'}};
par.RadPlot.p_vals       = Data.Behav.RT_TabSoloMEAN.pvalue';
par.RadPlot.maxValue     = maxRT;
par.RadPlot.minValue     = minRT;
RadPlot(RT_SoloS', RT_SoloK', par.RadPlot)

% ========================= PV =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'PV','Single');
par.RadPlot.unitMeasure  = 'pix/ms';
par.RadPlot.figureName   = SOLO_FIG;
par.RadPlot.cond_names1  = {'Act S'};
par.RadPlot.cond_names2  = {'Act K'};
par.RadPlot.InField      = 'PV';

par.RadPlot.color1       = Palette.color.ActS;
par.RadPlot.color2       = Palette.color.ActK;
par.RadPlot.lineStyle1   = LS_ACT;
par.RadPlot.lineStyle2   = LS_ACT;

par.RadPlot.legendNames  = {{'Monkey S','PV','Act S'}, ...
                            {'Monkey K','PV','Act K'}};
par.RadPlot.p_vals       = Data.Behav.PV_TabSoloMEAN.pvalue';
par.RadPlot.maxValue     = maxPV;
par.RadPlot.minValue     = minPV;
RadPlot(PV_SoloS', PV_SoloK', par.RadPlot)

% ========================= PVT =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'PVT','Single');
par.RadPlot.unitMeasure  = 'ms';
par.RadPlot.figureName   = SOLO_FIG;
par.RadPlot.cond_names1  = {'Act S'};
par.RadPlot.cond_names2  = {'Act K'};
par.RadPlot.InField      = 'PVT';

par.RadPlot.color1       = Palette.color.ActS;
par.RadPlot.color2       = Palette.color.ActK;
par.RadPlot.lineStyle1   = LS_ACT;
par.RadPlot.lineStyle2   = LS_ACT;

par.RadPlot.legendNames  = {{'Monkey S','PVT','Act S'}, ...
                            {'Monkey K','PVT','Act K'}};
par.RadPlot.p_vals       = Data.Behav.PVT_TabSoloMEAN.pvalue';
par.RadPlot.maxValue     = maxPVT;
par.RadPlot.minValue     = minPVT;
RadPlot(PVT_SoloS', PVT_SoloK', par.RadPlot)

%% =========================================================================
%% T-TEST: JOINT S vs JOINT K
%% =========================================================================

% ========================= AE =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'AE','Single');
par.RadPlot.unitMeasure  = '°';
par.RadPlot.figureName   = JOINT_FIG;
par.RadPlot.cond_names1  = {'Joint S'};
par.RadPlot.cond_names2  = {'Joint K'};
par.RadPlot.InField      = 'AE';

par.RadPlot.color1       = Palette.color.JointS;
par.RadPlot.color2       = Palette.color.JointK;
par.RadPlot.lineStyle1   = LS_JOINT;
par.RadPlot.lineStyle2   = LS_JOINT;

par.RadPlot.legendNames  = {{'Monkey S','AE','Joint S'}, ...
                            {'Monkey K','AE','Joint K'}};
par.RadPlot.p_vals       = Data.Behav.Jxy_TabJointMeanAE.pvalue';
par.RadPlot.maxValue     = maxAE;
par.RadPlot.minValue     = minAE;
RadPlot(AE_JointS', AE_JointK', par.RadPlot)

% ========================= AE_abs =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'AE_abs','Single');
par.RadPlot.unitMeasure  = '°';
par.RadPlot.figureName   = JOINT_FIG;
par.RadPlot.cond_names1  = {'Joint S'};
par.RadPlot.cond_names2  = {'Joint K'};
par.RadPlot.InField      = 'AE_abs';

par.RadPlot.color1       = Palette.color.JointS;
par.RadPlot.color2       = Palette.color.JointK;
par.RadPlot.lineStyle1   = LS_JOINT;
par.RadPlot.lineStyle2   = LS_JOINT;

par.RadPlot.legendNames  = {{'Monkey S','AE_abs','Joint S'}, ...
                            {'Monkey K','AE_abs','Joint K'}};
par.RadPlot.p_vals       = Data.Behav.Jxy_TabJointMeanAE_abs.pvalue';
par.RadPlot.maxValue     = maxAE_abs;
par.RadPlot.minValue     = minAE_abs;
RadPlot(AE_abs_JointS', AE_abs_JointK', par.RadPlot)

% ========================= ExT =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'ExT','Single');
par.RadPlot.unitMeasure  = 'ms';
par.RadPlot.figureName   = JOINT_FIG;
par.RadPlot.cond_names1  = {'Joint S'};
par.RadPlot.cond_names2  = {'Joint K'};
par.RadPlot.InField      = 'ExT';

par.RadPlot.color1       = Palette.color.JointS;
par.RadPlot.color2       = Palette.color.JointK;
par.RadPlot.lineStyle1   = LS_JOINT;
par.RadPlot.lineStyle2   = LS_JOINT;

par.RadPlot.legendNames  = {{'Monkey S','ExT','Joint S'}, ...
                            {'Monkey K','ExT','Joint K'}};
par.RadPlot.p_vals       = Data.Behav.Jxy_TabJointMeanExT.pvalue';
par.RadPlot.maxValue     = maxExT;
par.RadPlot.minValue     = minExT;
RadPlot(ExT_JointS', ExT_JointK', par.RadPlot)

% ========================= EC =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'EC','Single');
par.RadPlot.unitMeasure  = '%';
par.RadPlot.figureName   = JOINT_FIG;
par.RadPlot.cond_names1  = {'Joint S'};
par.RadPlot.cond_names2  = {'Joint K'};
par.RadPlot.InField      = 'EC';

par.RadPlot.color1       = Palette.color.JointS;
par.RadPlot.color2       = Palette.color.JointK;
par.RadPlot.lineStyle1   = LS_JOINT;
par.RadPlot.lineStyle2   = LS_JOINT;

par.RadPlot.legendNames  = {{'Monkey S','EC','Joint S'}, ...
                            {'Monkey K','EC','Joint K'}};
par.RadPlot.p_vals       = Data.Behav.Jxy_TabJointMeanEC.pvalue';
par.RadPlot.maxValue     = maxEC;
par.RadPlot.minValue     = minEC;
RadPlot(EC_JointS', EC_JointK', par.RadPlot)

% ========================= CMT =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'CMT','Single');
par.RadPlot.unitMeasure  = 'ms';
par.RadPlot.figureName   = JOINT_FIG;
par.RadPlot.cond_names1  = {'Joint S'};
par.RadPlot.cond_names2  = {'Joint K'};
par.RadPlot.InField      = 'CMT';

par.RadPlot.color1       = Palette.color.JointS;
par.RadPlot.color2       = Palette.color.JointK;
par.RadPlot.lineStyle1   = LS_JOINT;
par.RadPlot.lineStyle2   = LS_JOINT;

par.RadPlot.legendNames  = {{'Monkey S','CMT','Joint S'}, ...
                            {'Monkey K','CMT','Joint K'}};
par.RadPlot.p_vals       = Data.Behav.Jxy_TabJointMeanCMT.pvalue';
par.RadPlot.maxValue     = maxCMT;
par.RadPlot.minValue     = minCMT;
RadPlot(CMT_JointS', CMT_JointK', par.RadPlot)

% ========================= RT =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'RT','Single');
par.RadPlot.unitMeasure  = 'ms';
par.RadPlot.figureName   = JOINT_FIG;
par.RadPlot.cond_names1  = {'Joint S'};
par.RadPlot.cond_names2  = {'Joint K'};
par.RadPlot.InField      = 'RT';

par.RadPlot.color1       = Palette.color.JointS;
par.RadPlot.color2       = Palette.color.JointK;
par.RadPlot.lineStyle1   = LS_JOINT;
par.RadPlot.lineStyle2   = LS_JOINT;

par.RadPlot.legendNames  = {{'Monkey S','RT','Joint S'}, ...
                            {'Monkey K','RT','Joint K'}};
par.RadPlot.p_vals       = Data.Behav.RT_TabJointMEAN.pvalue';
par.RadPlot.maxValue     = maxRT;
par.RadPlot.minValue     = minRT;
RadPlot(RT_JointS', RT_JointK', par.RadPlot)

% ========================= PV =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'PV','Single');
par.RadPlot.unitMeasure  = 'pix/ms';
par.RadPlot.figureName   = JOINT_FIG;
par.RadPlot.cond_names1  = {'Joint S'};
par.RadPlot.cond_names2  = {'Joint K'};
par.RadPlot.InField      = 'PV';

par.RadPlot.color1       = Palette.color.JointS;
par.RadPlot.color2       = Palette.color.JointK;
par.RadPlot.lineStyle1   = LS_JOINT;
par.RadPlot.lineStyle2   = LS_JOINT;

par.RadPlot.legendNames  = {{'Monkey S','PV','Joint S'}, ...
                            {'Monkey K','PV','Joint K'}};
par.RadPlot.p_vals       = Data.Behav.PV_TabJointMEAN.pvalue';
par.RadPlot.maxValue     = maxPV;
par.RadPlot.minValue     = minPV;
RadPlot(PV_JointS', PV_JointK', par.RadPlot)

% ========================= PVT =========================
par.RadPlot.save_dir     = fullfile(ttestSaveDir,'PVT','Single');
par.RadPlot.unitMeasure  = 'ms';
par.RadPlot.figureName   = JOINT_FIG;
par.RadPlot.cond_names1  = {'Joint S'};
par.RadPlot.cond_names2  = {'Joint K'};
par.RadPlot.InField      = 'PVT';

par.RadPlot.color1       = Palette.color.JointS;
par.RadPlot.color2       = Palette.color.JointK;
par.RadPlot.lineStyle1   = LS_JOINT;
par.RadPlot.lineStyle2   = LS_JOINT;

par.RadPlot.legendNames  = {{'Monkey S','PVT','Joint S'}, ...
                            {'Monkey K','PVT','Joint K'}};
par.RadPlot.p_vals       = Data.Behav.PVT_TabJointMEAN.pvalue';
par.RadPlot.maxValue     = maxPVT;
par.RadPlot.minValue     = minPVT;
RadPlot(PVT_JointS', PVT_JointK', par.RadPlot)

%% =========================================================================
%% WITHIN-MONKEY RADAR PLOTS
%% =========================================================================

% ========================= AE =========================
par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'AE','Single');
par.RadPlotOwn.unitMeasure  = '°';
par.RadPlotOwn.InField       = 'AE';

par.RadPlotOwn.figureName    = 'ActS_JointS';
par.RadPlotOwn.cond_names1   = {'Act S'};
par.RadPlotOwn.cond_names2   = {'Joint S'};
par.RadPlotOwn.color1        = Palette.color.ActS;
par.RadPlotOwn.color2        = Palette.color.JointS;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;

par.RadPlotOwn.legendNames   = {{'Monkey S','AE','Act S'}, ...
                                {'Monkey S','AE','Joint S'}};
par.RadPlotOwn.p_vals        = Data.Behav.Jxy_TabMeanAE_S.pvalue';
par.RadPlotOwn.maxValue      = maxAE;
par.RadPlotOwn.minValue      = minAE;
RadPlotOwn(AE_SoloS', AE_JointS', par.RadPlotOwn)

par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'AE','Single');
par.RadPlotOwn.unitMeasure  = '°';
par.RadPlotOwn.InField       = 'AE';

par.RadPlotOwn.figureName    = 'ActK_JointK';
par.RadPlotOwn.cond_names1   = {'Act K'};
par.RadPlotOwn.cond_names2   = {'Joint K'};
par.RadPlotOwn.color1        = Palette.color.ActK;
par.RadPlotOwn.color2        = Palette.color.JointK;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;

par.RadPlotOwn.legendNames   = {{'Monkey K','AE','Act K'}, ...
                                {'Monkey K','AE','Joint K'}};
par.RadPlotOwn.p_vals        = Data.Behav.Jxy_TabMeanAE_K.pvalue';
par.RadPlotOwn.maxValue      = maxAE;
par.RadPlotOwn.minValue      = minAE;
RadPlotOwn(AE_SoloK', AE_JointK', par.RadPlotOwn)

% ========================= AE_abs =========================
par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'AE_abs','Single');
par.RadPlotOwn.unitMeasure  = '°';
par.RadPlotOwn.InField       = 'AE_abs';

par.RadPlotOwn.figureName    = 'ActS_JointS';
par.RadPlotOwn.cond_names1   = {'Act S'};
par.RadPlotOwn.cond_names2   = {'Joint S'};
par.RadPlotOwn.color1        = Palette.color.ActS;
par.RadPlotOwn.color2        = Palette.color.JointS;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;

par.RadPlotOwn.legendNames   = {{'Monkey S','AE_abs','Act S'}, ...
                                {'Monkey S','AE_abs','Joint S'}};
par.RadPlotOwn.p_vals        = Data.Behav.Jxy_TabMeanAE_abs_S.pvalue';
par.RadPlotOwn.maxValue      = maxAE_abs;
par.RadPlotOwn.minValue      = minAE_abs;
RadPlotOwn(AE_abs_SoloS', AE_abs_JointS', par.RadPlotOwn)

par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'AE_abs','Single');
par.RadPlotOwn.unitMeasure  = '°';
par.RadPlotOwn.InField       = 'AE_abs';

par.RadPlotOwn.figureName    = 'ActK_JointK';
par.RadPlotOwn.cond_names1   = {'Act K'};
par.RadPlotOwn.cond_names2   = {'Joint K'};
par.RadPlotOwn.color1        = Palette.color.ActK;
par.RadPlotOwn.color2        = Palette.color.JointK;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;

par.RadPlotOwn.legendNames   = {{'Monkey K','AE_abs','Act K'}, ...
                                {'Monkey K','AE_abs','Joint K'}};
par.RadPlotOwn.p_vals        = Data.Behav.Jxy_TabMeanAE_abs_K.pvalue';
par.RadPlotOwn.maxValue      = maxAE_abs;
par.RadPlotOwn.minValue      = minAE_abs;
RadPlotOwn(AE_abs_SoloK', AE_abs_JointK', par.RadPlotOwn)
% ========================= ExT =========================
par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'ExT','Single');
par.RadPlotOwn.unitMeasure   = 'ms';
par.RadPlotOwn.InField       = 'ExT';

par.RadPlotOwn.figureName    = 'ActS_JointS';
par.RadPlotOwn.cond_names1   = {'Act S'};
par.RadPlotOwn.cond_names2   = {'Joint S'};
par.RadPlotOwn.color1        = Palette.color.ActS;
par.RadPlotOwn.color2        = Palette.color.JointS;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;

par.RadPlotOwn.legendNames   = {{'Monkey S','ExT','Act S'}, ...
                                {'Monkey S','ExT','Joint S'}};
par.RadPlotOwn.p_vals        = Data.Behav.Jxy_TabMeanExT_S.pvalue';
par.RadPlotOwn.maxValue      = maxExT;
par.RadPlotOwn.minValue      = minExT;
RadPlotOwn(ExT_SoloS', ExT_JointS', par.RadPlotOwn)

par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'ExT','Single');
par.RadPlotOwn.unitMeasure   = 'ms';
par.RadPlotOwn.InField       = 'ExT';

par.RadPlotOwn.figureName    = 'ActK_JointK';
par.RadPlotOwn.cond_names1   = {'Act K'};
par.RadPlotOwn.cond_names2   = {'Joint K'};
par.RadPlotOwn.color1        = Palette.color.ActK;
par.RadPlotOwn.color2        = Palette.color.JointK;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;

par.RadPlotOwn.legendNames   = {{'Monkey K','ExT','Act K'}, ...
                                {'Monkey K','ExT','Joint K'}};
par.RadPlotOwn.p_vals        = Data.Behav.Jxy_TabMeanExT_K.pvalue';
par.RadPlotOwn.maxValue      = maxExT;
par.RadPlotOwn.minValue      = minExT;
RadPlotOwn(ExT_SoloK', ExT_JointK', par.RadPlotOwn)

% ========================= EC =========================
par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'EC','Single');
par.RadPlotOwn.unitMeasure   = '%';
par.RadPlotOwn.InField       = 'EC';

par.RadPlotOwn.figureName    = 'ActS_JointS';
par.RadPlotOwn.cond_names1   = {'Act S'};
par.RadPlotOwn.cond_names2   = {'Joint S'};
par.RadPlotOwn.color1        = Palette.color.ActS;
par.RadPlotOwn.color2        = Palette.color.JointS;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;

par.RadPlotOwn.legendNames   = {{'Monkey S','EC','Act S'}, ...
                                {'Monkey S','EC','Joint S'}};
par.RadPlotOwn.p_vals        = Data.Behav.Jxy_TabMeanEC_S.pvalue';
par.RadPlotOwn.maxValue      = maxEC;
par.RadPlotOwn.minValue      = minEC;
RadPlotOwn(EC_SoloS', EC_JointS', par.RadPlotOwn)

par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'EC','Single');
par.RadPlotOwn.unitMeasure   = '%';
par.RadPlotOwn.InField       = 'EC';

par.RadPlotOwn.figureName    = 'ActK_JointK';
par.RadPlotOwn.cond_names1   = {'Act K'};
par.RadPlotOwn.cond_names2   = {'Joint K'};
par.RadPlotOwn.color1        = Palette.color.ActK;
par.RadPlotOwn.color2        = Palette.color.JointK;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;

par.RadPlotOwn.legendNames   = {{'Monkey K','EC','Act K'}, ...
                                {'Monkey K','EC','Joint K'}};
par.RadPlotOwn.p_vals        = Data.Behav.Jxy_TabMeanEC_K.pvalue';
par.RadPlotOwn.maxValue      = maxEC;
par.RadPlotOwn.minValue      = minEC;
RadPlotOwn(EC_SoloK', EC_JointK', par.RadPlotOwn)

% ========================= CMT =========================
par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'CMT','Single');
par.RadPlotOwn.unitMeasure   = 'ms';
par.RadPlotOwn.InField       = 'CMT';

par.RadPlotOwn.figureName    = 'ActS_JointS';
par.RadPlotOwn.cond_names1   = {'Act S'};
par.RadPlotOwn.cond_names2   = {'Joint S'};
par.RadPlotOwn.color1        = Palette.color.ActS;
par.RadPlotOwn.color2        = Palette.color.JointS;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;
par.RadPlotOwn.legendNames   = {{'Monkey S','CMT','Act S'}, ...
                                {'Monkey S','CMT','Joint S'}};
par.RadPlotOwn.p_vals        = Data.Behav.Jxy_TabMeanCMT_S.pvalue';
par.RadPlotOwn.maxValue      = maxCMT;
par.RadPlotOwn.minValue      = minCMT;
RadPlotOwn(CMT_SoloS', CMT_JointS', par.RadPlotOwn)
par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'CMT','Single');
par.RadPlotOwn.unitMeasure   = 'ms';
par.RadPlotOwn.InField       = 'CMT';
par.RadPlotOwn.figureName    = 'ActK_JointK';
par.RadPlotOwn.cond_names1   = {'Act K'};
par.RadPlotOwn.cond_names2   = {'Joint K'};
par.RadPlotOwn.color1        = Palette.color.ActK;
par.RadPlotOwn.color2        = Palette.color.JointK;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;
par.RadPlotOwn.legendNames   = {{'Monkey K','CMT','Act K'}, ...
                                {'Monkey K','CMT','Joint K'}};
par.RadPlotOwn.p_vals        = Data.Behav.Jxy_TabMeanCMT_K.pvalue';
par.RadPlotOwn.maxValue      = maxCMT;
par.RadPlotOwn.minValue      = minCMT;
RadPlotOwn(CMT_SoloK', CMT_JointK', par.RadPlotOwn)

% ========================= RT =========================
par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'RT','Single');
par.RadPlotOwn.unitMeasure   = 'ms';
par.RadPlotOwn.InField       = 'RT';

par.RadPlotOwn.figureName    = 'ActS_JointS';
par.RadPlotOwn.cond_names1   = {'Act S'};
par.RadPlotOwn.cond_names2   = {'Joint S'};
par.RadPlotOwn.color1        = Palette.color.ActS;
par.RadPlotOwn.color2        = Palette.color.JointS;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;
par.RadPlotOwn.legendNames   = {{'Monkey S','RT','Act S'}, ...
                                {'Monkey S','RT','Joint S'}};
par.RadPlotOwn.p_vals        = Data.Behav.RT_TabMEAN_S.pvalue';
par.RadPlotOwn.maxValue      = maxRT;
par.RadPlotOwn.minValue      = minRT;
RadPlotOwn(RT_SoloS', RT_JointS', par.RadPlotOwn)

par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'RT','Single');
par.RadPlotOwn.unitMeasure   = 'ms';
par.RadPlotOwn.InField       = 'RT';
par.RadPlotOwn.figureName    = 'ActK_JointK';
par.RadPlotOwn.cond_names1   = {'Act K'};
par.RadPlotOwn.cond_names2   = {'Joint K'};
par.RadPlotOwn.color1        = Palette.color.ActK;
par.RadPlotOwn.color2        = Palette.color.JointK;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;
par.RadPlotOwn.legendNames   = {{'Monkey K','RT','Act K'}, ...
                                {'Monkey K','RT','Joint K'}};
par.RadPlotOwn.p_vals        = Data.Behav.RT_TabMEAN_K.pvalue';
par.RadPlotOwn.maxValue      = maxRT;
par.RadPlotOwn.minValue      = minRT;
RadPlotOwn(RT_SoloK', RT_JointK', par.RadPlotOwn)

% ========================= PV =========================
par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'PV','Single');
par.RadPlotOwn.unitMeasure   = 'pix/ms';
par.RadPlotOwn.InField       = 'PV';
par.RadPlotOwn.figureName    = 'ActS_JointS';
par.RadPlotOwn.cond_names1   = {'Act S'};
par.RadPlotOwn.cond_names2   = {'Joint S'};
par.RadPlotOwn.color1        = Palette.color.ActS;
par.RadPlotOwn.color2        = Palette.color.JointS;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;
par.RadPlotOwn.legendNames   = {{'Monkey S','PV','Act S'}, ...
                                {'Monkey S','PV','Joint S'}};
par.RadPlotOwn.p_vals        = Data.Behav.PV_TabMEAN_S.pvalue';
par.RadPlotOwn.maxValue      = maxPV;
par.RadPlotOwn.minValue      = minPV;
RadPlotOwn(PV_SoloS', PV_JointS', par.RadPlotOwn)

par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'PV','Single');
par.RadPlotOwn.unitMeasure   = 'pix/ms';
par.RadPlotOwn.InField       = 'PV';
par.RadPlotOwn.figureName    = 'ActK_JointK';
par.RadPlotOwn.cond_names1   = {'Act K'};
par.RadPlotOwn.cond_names2   = {'Joint K'};
par.RadPlotOwn.color1        = Palette.color.ActK;
par.RadPlotOwn.color2        = Palette.color.JointK;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;
par.RadPlotOwn.legendNames   = {{'Monkey K','PV','Act K'}, ...
                                {'Monkey K','PV','Joint K'}};
par.RadPlotOwn.p_vals        = Data.Behav.PV_TabMEAN_K.pvalue';
par.RadPlotOwn.maxValue      = maxPV;
par.RadPlotOwn.minValue      = minPV;
RadPlotOwn(PV_SoloK', PV_JointK', par.RadPlotOwn)

% ========================= PVT =========================
par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'PVT','Single');
par.RadPlotOwn.unitMeasure   = 'ms';
par.RadPlotOwn.InField       = 'PVT';
par.RadPlotOwn.figureName    = 'ActS_JointS';
par.RadPlotOwn.cond_names1   = {'Act S'};
par.RadPlotOwn.cond_names2   = {'Joint S'};
par.RadPlotOwn.color1        = Palette.color.ActS;
par.RadPlotOwn.color2        = Palette.color.JointS;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;
par.RadPlotOwn.legendNames   = {{'Monkey S','PVT','Act S'}, ...
                                {'Monkey S','PVT','Joint S'}};
par.RadPlotOwn.p_vals        = Data.Behav.PVT_TabMEAN_S.pvalue';
par.RadPlotOwn.maxValue      = maxPVT;
par.RadPlotOwn.minValue      = minPVT;
RadPlotOwn(PVT_SoloS', PVT_JointS', par.RadPlotOwn)

par.RadPlotOwn.save_dir      = fullfile(ttestSaveDir,'PVT','Single');
par.RadPlotOwn.unitMeasure   = 'ms';
par.RadPlotOwn.InField       = 'PVT';
par.RadPlotOwn.figureName    = 'ActK_JointK';
par.RadPlotOwn.cond_names1   = {'Act K'};
par.RadPlotOwn.cond_names2   = {'Joint K'};
par.RadPlotOwn.color1        = Palette.color.ActK;
par.RadPlotOwn.color2        = Palette.color.JointK;
par.RadPlotOwn.lineStyle1    = LS_ACT;
par.RadPlotOwn.lineStyle2    = LS_JOINT;
par.RadPlotOwn.legendNames   = {{'Monkey K','PVT','Act K'}, ...
                                {'Monkey K','PVT','Joint K'}};
par.RadPlotOwn.p_vals        = Data.Behav.PVT_TabMEAN_K.pvalue';
par.RadPlotOwn.maxValue      = maxPVT;
par.RadPlotOwn.minValue      = minPVT;
RadPlotOwn(PVT_SoloK', PVT_JointK', par.RadPlotOwn)

%% =========================================================================
%% FOUR PLOT: ACT/JOINT S-K
%% =========================================================================

% ========================= AE =========================
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'AE','Multiple');
par.RadMultiplePlot.InField       = 'AE';
par.RadMultiplePlot.unitMeasure   = '°';
par.RadMultiplePlot.figureName    = 'ActS_JointS_ActK_JointK';
par.RadMultiplePlot.cond_names1   = {'Act S'};
par.RadMultiplePlot.cond_names2   = {'Joint S'};
par.RadMultiplePlot.cond_names3   = {'Act K'};
par.RadMultiplePlot.cond_names4   = {'Joint K'};
par.RadMultiplePlot.color1        = Palette.color.ActS;
par.RadMultiplePlot.color2        = Palette.color.JointS;
par.RadMultiplePlot.color3        = Palette.color.ActK;
par.RadMultiplePlot.color4        = Palette.color.JointK;
par.RadMultiplePlot.lineStyle1    = LS_ACT;
par.RadMultiplePlot.lineStyle2    = LS_JOINT;
par.RadMultiplePlot.lineStyle3    = LS_ACT;
par.RadMultiplePlot.lineStyle4    = LS_JOINT;
par.RadMultiplePlot.legendNames   = {...
    {'Monkey S','AE','Act S'}, ...
    {'Monkey S','AE','Joint S'}, ...
    {'Monkey K','AE','Act K'}, ...
    {'Monkey K','AE','Joint K'}};
par.RadMultiplePlot.maxValue      = maxAE;
par.RadMultiplePlot.minValue      = minAE;
RadMultiplePlot(AE_SoloS', AE_JointS', AE_SoloK', AE_JointK', par.RadMultiplePlot);

% ========================= AE_abs =========================
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'AE_abs','Multiple');
par.RadMultiplePlot.InField       = 'AE_abs';
par.RadMultiplePlot.unitMeasure   = '°';
par.RadMultiplePlot.figureName    = 'ActS_JointS_ActK_JointK';
par.RadMultiplePlot.cond_names1   = {'Act S'};
par.RadMultiplePlot.cond_names2   = {'Joint S'};
par.RadMultiplePlot.cond_names3   = {'Act K'};
par.RadMultiplePlot.cond_names4   = {'Joint K'};
par.RadMultiplePlot.color1        = Palette.color.ActS;
par.RadMultiplePlot.color2        = Palette.color.JointS;
par.RadMultiplePlot.color3        = Palette.color.ActK;
par.RadMultiplePlot.color4        = Palette.color.JointK;
par.RadMultiplePlot.lineStyle1    = LS_ACT;
par.RadMultiplePlot.lineStyle2    = LS_JOINT;
par.RadMultiplePlot.lineStyle3    = LS_ACT;
par.RadMultiplePlot.lineStyle4    = LS_JOINT;
par.RadMultiplePlot.legendNames   = {...
    {'Monkey S','AE_abs','Act S'}, ...
    {'Monkey S','AE_abs','Joint S'}, ...
    {'Monkey K','AE_abs','Act K'}, ...
    {'Monkey K','AE_abs','Joint K'}};
par.RadMultiplePlot.maxValue      = maxAE_abs;
par.RadMultiplePlot.minValue      = minAE_abs;
RadMultiplePlot(AE_abs_SoloS', AE_abs_JointS', AE_abs_SoloK', AE_abs_JointK', par.RadMultiplePlot);

% ========================= ExT =========================
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'ExT','Multiple');
par.RadMultiplePlot.unitMeasure   = 'ms';
par.RadMultiplePlot.InField       = 'ExT';
par.RadMultiplePlot.figureName    = 'ActS_JointS_ActK_JointK';
par.RadMultiplePlot.cond_names1   = {'Act S'};
par.RadMultiplePlot.cond_names2   = {'Joint S'};
par.RadMultiplePlot.cond_names3   = {'Act K'};
par.RadMultiplePlot.cond_names4   = {'Joint K'};
par.RadMultiplePlot.color1        = Palette.color.ActS;
par.RadMultiplePlot.color2        = Palette.color.JointS;
par.RadMultiplePlot.color3        = Palette.color.ActK;
par.RadMultiplePlot.color4        = Palette.color.JointK;
par.RadMultiplePlot.lineStyle1    = LS_ACT;
par.RadMultiplePlot.lineStyle2    = LS_JOINT;
par.RadMultiplePlot.lineStyle3    = LS_ACT;
par.RadMultiplePlot.lineStyle4    = LS_JOINT;
par.RadMultiplePlot.legendNames   = {...
    {'Monkey S','ExT','Act S'}, ...
    {'Monkey S','ExT','Joint S'}, ...
    {'Monkey K','ExT','Act K'}, ...
    {'Monkey K','ExT','Joint K'}};
par.RadMultiplePlot.maxValue      = maxExT;
par.RadMultiplePlot.minValue      = minExT;
RadMultiplePlot(ExT_SoloS', ExT_JointS', ExT_SoloK', ExT_JointK', par.RadMultiplePlot);

% ========================= EC =========================
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'EC','Multiple');
par.RadMultiplePlot.unitMeasure   = '%';
par.RadMultiplePlot.InField       = 'EC';
par.RadMultiplePlot.figureName    = 'ActS_JointS_ActK_JointK';
par.RadMultiplePlot.cond_names1   = {'Act S'};
par.RadMultiplePlot.cond_names2   = {'Joint S'};
par.RadMultiplePlot.cond_names3   = {'Act K'};
par.RadMultiplePlot.cond_names4   = {'Joint K'};
par.RadMultiplePlot.color1        = Palette.color.ActS;
par.RadMultiplePlot.color2        = Palette.color.JointS;
par.RadMultiplePlot.color3        = Palette.color.ActK;
par.RadMultiplePlot.color4        = Palette.color.JointK;
par.RadMultiplePlot.lineStyle1    = LS_ACT;
par.RadMultiplePlot.lineStyle2    = LS_JOINT;
par.RadMultiplePlot.lineStyle3    = LS_ACT;
par.RadMultiplePlot.lineStyle4    = LS_JOINT;
par.RadMultiplePlot.legendNames   = {...
    {'Monkey S','EC','Act S'}, ...
    {'Monkey S','EC','Joint S'}, ...
    {'Monkey K','EC','Act K'}, ...
    {'Monkey K','EC','Joint K'}};
par.RadMultiplePlot.maxValue      = maxEC;
par.RadMultiplePlot.minValue      = minEC;
RadMultiplePlot(EC_SoloS', EC_JointS', EC_SoloK', EC_JointK', par.RadMultiplePlot);

% ========================= CMT =========================
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'CMT','Multiple');
par.RadMultiplePlot.unitMeasure   = 'ms';
par.RadMultiplePlot.InField       = 'CMT';

par.RadMultiplePlot.figureName    = 'ActS_JointS_ActK_JointK';
par.RadMultiplePlot.cond_names1   = {'Act S'};
par.RadMultiplePlot.cond_names2   = {'Joint S'};
par.RadMultiplePlot.cond_names3   = {'Act K'};
par.RadMultiplePlot.cond_names4   = {'Joint K'};
par.RadMultiplePlot.color1        = Palette.color.ActS;
par.RadMultiplePlot.color2        = Palette.color.JointS;
par.RadMultiplePlot.color3        = Palette.color.ActK;
par.RadMultiplePlot.color4        = Palette.color.JointK;
par.RadMultiplePlot.lineStyle1    = LS_ACT;
par.RadMultiplePlot.lineStyle2    = LS_JOINT;
par.RadMultiplePlot.lineStyle3    = LS_ACT;
par.RadMultiplePlot.lineStyle4    = LS_JOINT;
par.RadMultiplePlot.legendNames   = {...
    {'Monkey S','CMT','Act S'}, ...
    {'Monkey S','CMT','Joint S'}, ...
    {'Monkey K','CMT','Act K'}, ...
    {'Monkey K','CMT','Joint K'}};
par.RadMultiplePlot.maxValue      = maxCMT;
par.RadMultiplePlot.minValue      = minCMT;
RadMultiplePlot(CMT_SoloS', CMT_JointS', CMT_SoloK', CMT_JointK', par.RadMultiplePlot);

% ========================= RT =========================
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'RT','Multiple');
par.RadMultiplePlot.unitMeasure   = 'ms';
par.RadMultiplePlot.InField       = 'RT';

par.RadMultiplePlot.figureName    = 'ActS_JointS_ActK_JointK';
par.RadMultiplePlot.cond_names1   = {'Act S'};
par.RadMultiplePlot.cond_names2   = {'Joint S'};
par.RadMultiplePlot.cond_names3   = {'Act K'};
par.RadMultiplePlot.cond_names4   = {'Joint K'};
par.RadMultiplePlot.color1        = Palette.color.ActS;
par.RadMultiplePlot.color2        = Palette.color.JointS;
par.RadMultiplePlot.color3        = Palette.color.ActK;
par.RadMultiplePlot.color4        = Palette.color.JointK;
par.RadMultiplePlot.lineStyle1    = LS_ACT;
par.RadMultiplePlot.lineStyle2    = LS_JOINT;
par.RadMultiplePlot.lineStyle3    = LS_ACT;
par.RadMultiplePlot.lineStyle4    = LS_JOINT;

par.RadMultiplePlot.legendNames   = {...
    {'Monkey S','RT','Act S'}, ...
    {'Monkey S','RT','Joint S'}, ...
    {'Monkey K','RT','Act K'}, ...
    {'Monkey K','RT','Joint K'}};
par.RadMultiplePlot.maxValue      = maxRT;
par.RadMultiplePlot.minValue      = minRT;
RadMultiplePlot(RT_SoloS', RT_JointS', RT_SoloK', RT_JointK', par.RadMultiplePlot);

% ========================= PV =========================
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'PV','Multiple');
par.RadMultiplePlot.unitMeasure   = 'pix/ms';
par.RadMultiplePlot.InField       = 'PV';

par.RadMultiplePlot.figureName    = 'ActS_JointS_ActK_JointK';
par.RadMultiplePlot.cond_names1   = {'Act S'};
par.RadMultiplePlot.cond_names2   = {'Joint S'};
par.RadMultiplePlot.cond_names3   = {'Act K'};
par.RadMultiplePlot.cond_names4   = {'Joint K'};
par.RadMultiplePlot.color1        = Palette.color.ActS;
par.RadMultiplePlot.color2        = Palette.color.JointS;
par.RadMultiplePlot.color3        = Palette.color.ActK;
par.RadMultiplePlot.color4        = Palette.color.JointK;
par.RadMultiplePlot.lineStyle1    = LS_ACT;
par.RadMultiplePlot.lineStyle2    = LS_JOINT;
par.RadMultiplePlot.lineStyle3    = LS_ACT;
par.RadMultiplePlot.lineStyle4    = LS_JOINT;

par.RadMultiplePlot.legendNames   = {...
    {'Monkey S','PV','Act S'},...
    {'Monkey S','PV','Joint S'},...
    {'Monkey K','PV','Act K'},...
    {'Monkey K','PV','Joint K'}};
par.RadMultiplePlot.maxValue      = maxPV;
par.RadMultiplePlot.minValue      = minPV;
RadMultiplePlot(PV_SoloS', PV_JointS', PV_SoloK',PV_JointK', par.RadMultiplePlot);

% ========================= PVT =========================
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'PVT','Multiple');
par.RadMultiplePlot.unitMeasure   = 'ms';
par.RadMultiplePlot.InField       = 'PVT';

par.RadMultiplePlot.figureName    = 'ActS_JointS_ActK_JointK';
par.RadMultiplePlot.cond_names1   = {'Act S'};
par.RadMultiplePlot.cond_names2   = {'Joint S'};
par.RadMultiplePlot.cond_names3   = {'Act K'};
par.RadMultiplePlot.cond_names4   = {'Joint K'};
par.RadMultiplePlot.color1        = Palette.color.ActS;
par.RadMultiplePlot.color2        = Palette.color.JointS;
par.RadMultiplePlot.color3        = Palette.color.ActK;
par.RadMultiplePlot.color4        = Palette.color.JointK;
par.RadMultiplePlot.lineStyle1    = LS_ACT;
par.RadMultiplePlot.lineStyle2    = LS_JOINT;
par.RadMultiplePlot.lineStyle3    = LS_ACT;
par.RadMultiplePlot.lineStyle4    = LS_JOINT;

par.RadMultiplePlot.legendNames   = {...
    {'Monkey S','PVT','Act S'},...
    {'Monkey S','PVT','Joint S'},...
    {'Monkey K','PVT','Act K'},...
    {'Monkey K','PVT','Joint K'}};
par.RadMultiplePlot.maxValue      = maxPVT;
par.RadMultiplePlot.minValue      = minPVT;
RadMultiplePlot(PVT_SoloS', PVT_JointS', PVT_SoloK',PVT_JointK', par.RadMultiplePlot);

%% =========================================================================
%% ICD
%% =========================================================================
ICDmax  = icd.Data.Behav.ICDmax_total_Mean;
ICDmean = icd.Data.Behav.ICDmean_total_Mean;
ICDauc  = icd.Data.Behav.ICDauc_total_Mean;

ICDmax_Joint   = ICDmax(3,:);
ICDmean_Joint  = ICDmean(3,:);
ICDauc_Joint   = ICDauc(3,:);

frac = 0.15;

deltaICDmax    = frac * max(ICDmax_Joint, [], 'omitnan');
maxICDmax      = max(ICDmax_Joint, [], 'omitnan') + deltaICDmax;
minICDmax      = min(ICDmax_Joint, [], 'omitnan') - deltaICDmax;

deltaICDmean   = frac * max(ICDmean_Joint, [], 'omitnan');
maxICDmean     = max(ICDmean_Joint, [], 'omitnan') + deltaICDmean;
minICDmean     = min(ICDmean_Joint, [], 'omitnan') - deltaICDmean;

deltaICDauc    = frac * max(ICDauc_Joint, [], 'omitnan');
maxICDauc      = max(ICDauc_Joint, [], 'omitnan') + deltaICDauc;
minICDauc      = min(ICDauc_Joint, [], 'omitnan') - deltaICDauc;

% ICD max
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'ICD');
par.RadMultiplePlot.unitMeasure   = 'dva';
par.RadMultiplePlot.InField       = 'ICDmax';

par.RadMultiplePlot.figureName    = 'Joint';
par.RadMultiplePlot.cond_names1   = {'ICD max'};
par.RadMultiplePlot.cond_names2   = {''};
par.RadMultiplePlot.cond_names3   = {''};
par.RadMultiplePlot.cond_names4   = {''};
par.RadMultiplePlot.color1        = [0 0 0]; %[0.22, 0.49, 0.72];
par.RadMultiplePlot.color2        = [];
par.RadMultiplePlot.color3        = [];
par.RadMultiplePlot.color4        = [];
par.RadMultiplePlot.lineStyle1    = '-.';
par.RadMultiplePlot.lineStyle2    = '';
par.RadMultiplePlot.lineStyle3    = '';
par.RadMultiplePlot.lineStyle4    = '';

par.RadMultiplePlot.legendNames   = {{'Monkey S and K','ICD max','Joint'}};
par.RadMultiplePlot.maxValue      = maxICDmax;
par.RadMultiplePlot.minValue      = minICDmax;
RadMultiplePlot(ICDmax_Joint, par.RadMultiplePlot);

% ICD mean
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'ICD');
par.RadMultiplePlot.unitMeasure   = 'dva';
par.RadMultiplePlot.InField       = 'ICDmean';

par.RadMultiplePlot.figureName    = 'Joint';
par.RadMultiplePlot.cond_names1   = {'ICD mean'};
par.RadMultiplePlot.cond_names2   = {''};
par.RadMultiplePlot.cond_names3   = {''};
par.RadMultiplePlot.cond_names4   = {''};
par.RadMultiplePlot.color1        =  [0 0 0]; %[0.60, 0.31, 0.64];
par.RadMultiplePlot.color2        = [];
par.RadMultiplePlot.color3        = [];
par.RadMultiplePlot.color4        = [];
par.RadMultiplePlot.lineStyle1    = '-.';
par.RadMultiplePlot.lineStyle2    = '';
par.RadMultiplePlot.lineStyle3    = '';
par.RadMultiplePlot.lineStyle4    = '';

par.RadMultiplePlot.legendNames   = {{'Monkey S and K','ICD mean','Joint'}};
par.RadMultiplePlot.maxValue      = maxICDmean;
par.RadMultiplePlot.minValue      = minICDmean;
RadMultiplePlot(ICDmean_Joint, par.RadMultiplePlot);

% ICD auc
par.RadMultiplePlot.save_dir      = fullfile(ttestSaveDir,'ICD');
par.RadMultiplePlot.unitMeasure   = 'dva';
par.RadMultiplePlot.InField       = 'ICDauc';

par.RadMultiplePlot.figureName    = 'Joint';
par.RadMultiplePlot.cond_names1   = {'ICD auc'};
par.RadMultiplePlot.cond_names2   = {''};
par.RadMultiplePlot.cond_names3   = {''};
par.RadMultiplePlot.cond_names4   = {''};
par.RadMultiplePlot.color1        = [0 0 0]; %[0.00, 0.62, 0.45];
par.RadMultiplePlot.color2        = [];
par.RadMultiplePlot.color3        = [];
par.RadMultiplePlot.color4        = [];
par.RadMultiplePlot.lineStyle1    = '-.';
par.RadMultiplePlot.lineStyle2    = '';
par.RadMultiplePlot.lineStyle3    = '';
par.RadMultiplePlot.lineStyle4    = '';

par.RadMultiplePlot.legendNames   = {{'Monkey S and K','ICD auc','Joint'}};
par.RadMultiplePlot.maxValue      = maxICDauc;
par.RadMultiplePlot.minValue      = minICDauc;
RadMultiplePlot(ICDauc_Joint, par.RadMultiplePlot);

set(0,'DefaultFigureVisible',origState);