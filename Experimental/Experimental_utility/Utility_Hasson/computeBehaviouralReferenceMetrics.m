function behavMetrics = computeBehaviouralReferenceMetrics(Data)
% COMPUTEBEHAVIOURALREFERENCEMETRICS
%
% Purpose
%   Computes behavioural reference metrics used for the interpretation
%   and annotation of Hasson directional coupling analyses.
%
%   Reaction time (RT) and exit time (ExT) measurements are extracted
%   from the behavioural dataset and averaged across monkeys and task
%   conditions to obtain global temporal reference markers.
%
% Inputs
%   Data
%       Behavioural structure loaded from the session summary file.
%
% Outputs
%   behavMetrics
%       Structure containing:
%
%       .RT_mean
%           Mean reaction time across monkeys and conditions [ms]
%
%       .ExT_mean
%           Mean exit time across monkeys and conditions [ms]
%
% Used by
%   MAIN_HASSON_WITH_EXTRACT_SPIKES
%
% Author
%   Mirco Frosolone
%

Jxy_S    = Data.Behav.Jxy_S;
Jxy_K    = Data.Behav.Jxy_K;
RT_total = Data.Behav.RT_total;
behavMetrics = struct();
%% Exit times

[~,~,ExT] = extractBehavDeltas( ...
    Jxy_S,...
    Jxy_K,...
    'exit_time');

ExT_SoloK  = mean(ExT.K.Solo ,'omitnan')';
ExT_SoloS  = mean(ExT.S.Solo ,'omitnan')';
ExT_JointK = mean(ExT.K.Joint,'omitnan')';
ExT_JointS = mean(ExT.S.Joint,'omitnan')';

ExT_mean = 1000 * mean( ...
    [ExT_SoloK;
     ExT_SoloS;
     ExT_JointK;
     ExT_JointS], ...
    'omitnan');

%% Reaction times

[~,~,RT] = extractFromTotalStruct( ...
    RT_total,...
    'less');

RT_SoloK  = mean(RT.K.Solo ,'omitnan')';
RT_SoloS  = mean(RT.S.Solo ,'omitnan')';
RT_JointK = mean(RT.K.Joint,'omitnan')';
RT_JointS = mean(RT.S.Joint,'omitnan')';

RT_mean = mean( ...
    [RT_SoloK;
     RT_SoloS;
     RT_JointK;
     RT_JointS], ...
    'omitnan');
behavMetrics.RT_mean = RT_mean;
behavMetrics.ExT_mean = ExT_mean;
end