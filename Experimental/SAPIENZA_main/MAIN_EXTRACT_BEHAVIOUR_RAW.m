%% MAIN_EXTRACT_BEHAVIOUR_RAW
% Extract behavioural metrics from RAW sessions and save them by chamber/session.
%
% Output layout:
%   /TESTS/SAPIENZA/BEHAVIOUR/
%       Frontal/SK007/Behaviour_RAW_FrontalSK007_D1D2D3D4D5D6D7D8_C1C2C3.mat
%       Frontal/Behaviour_RAW_Frontal_AllSessions_D1D2D3D4D5D6D7D8_C1C2C3.mat
%       Parietal/...
%       Behaviour_RAW_AllChambers_D1D2D3D4D5D6D7D8_C1C2C3.mat

currentScriptDir = fileparts(mfilename('fullpath'));
trialBoxRoot = fileparts(fileparts(currentScriptDir));
addpath(fullfile(trialBoxRoot, 'Experimental', 'Experimental_utility', 'Utility_Behaviour'));
addpath(fullfile(trialBoxRoot, 'Experimental', 'SAPIENZA_main'));

outputRoot = '/TESTS/SAPIENZA/BEHAVIOUR';
fileConditionTag = 'D1D2D3D4D5D6D7D8_C1C2C3';

extractionPms = struct();
extractionPms.selectionLFPs = 'maximum';

%% LFP SESSIONS
% FRONTAL
sessionNamesFrontal = {...
    'SK007', 'SK008', 'SK009', 'SK011', 'SK012', 'SK013', 'SK020', 'SK021', ...
    'SK022', 'SK023', 'SK024', 'SK025', 'SK026', 'SK028', 'SK029', ...
    'SK030', 'SK031', 'SK032', 'SK033', 'SK035', 'SK036', 'SK037', 'SK038', ...
    'SK039', 'SK040', 'SK041', 'SK042', 'SK043', 'SK044', 'SK045', 'SK046', ...
    'SK047', 'SK048'};

% PARIETAL
sessionNamesParietal = {...
    'SK004', 'SK005', 'SK049', 'SK050', 'SK051', 'SK052', 'SK054', 'SK055', ...
    'SK056', 'SK057', 'SK059', 'SK060', 'SK061', 'SK062', 'SK063', 'SK064', ...
    'SK065', 'SK066', 'SK067', 'SK069', 'SK070', 'SK071', 'SK072', 'SK074', ...
    'SK075'};

sessionNamesByChamber = struct();
sessionNamesByChamber.Frontal = sessionNamesFrontal;
sessionNamesByChamber.Parietal = sessionNamesParietal;
chamberNames = fieldnames(sessionNamesByChamber);

conditionIndexMap = buildBehaviourConditionIndexMap();

behaviourByChamber = struct();
metadataByChamber = struct();
summaryRows = cell(0, 5);

if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

for chamberCounter = 1:numel(chamberNames)
    chamberName = chamberNames{chamberCounter};
    sessionNames = sessionNamesByChamber.(chamberName);
    chamberOutputDir = fullfile(outputRoot, chamberName);

    if ~exist(chamberOutputDir, 'dir')
        mkdir(chamberOutputDir);
    end

    chamberBehaviour = struct();
    chamberMetadata = struct();
    chamberMetadata.chamber = chamberName;
    chamberMetadata.sessionNames = sessionNames;
    chamberMetadata.outputDir = chamberOutputDir;
    chamberMetadata.extractionPms = extractionPms;
    chamberMetadata.conditionIndexMap = conditionIndexMap;
    chamberMetadata.conditionIndexDescription = [
        "1:8 = Act S / Obs K, directions 1:8"; ...
        "9:16 = Obs S / Act K, directions 1:8"; ...
        "17:24 = Act S / Act K, directions 1:8"];
    chamberMetadata.failedSessions = struct('sessionName', {}, 'message', {});

    fprintf('Extracting RAW behaviour for %s sessions...\n', chamberName);

    for sessionCounter = 1:numel(sessionNames)
        sessionName = sessionNames{sessionCounter};
        sessionOutputDir = fullfile(chamberOutputDir, sessionName);

        if ~exist(sessionOutputDir, 'dir')
            mkdir(sessionOutputDir);
        end

        sessionFileName = sprintf('Behaviour_RAW_%s%s_%s.mat', chamberName, sessionName, fileConditionTag);
        sessionOutputPath = fullfile(sessionOutputDir, sessionFileName);

        fprintf('  %s %s (%d/%d)\n', chamberName, sessionName, sessionCounter, numel(sessionNames));

        try
            behaviourData = extractBehaviourLFP(sessionName, extractionPms);
            sessionMetadata = chamberMetadata;
            sessionMetadata.sessionName = sessionName;
            sessionMetadata.sessionOutputDir = sessionOutputDir;
            sessionMetadata.sessionOutputPath = sessionOutputPath;

            Data = struct();
            Data.Behav = behaviourData;
            Data.meta = sessionMetadata;

            save(sessionOutputPath, 'Data', 'behaviourData', 'sessionMetadata', '-v7.3');
            chamberBehaviour.(sessionName) = behaviourData;
            summaryRows(end + 1, :) = {chamberName, sessionName, true, sessionOutputPath, ''}; %#ok<SAGROW>
        catch extractionError
            warning('MAIN_EXTRACT_BEHAVIOUR_RAW:SessionFailed', ...
                'Skipping %s %s: %s', chamberName, sessionName, extractionError.message);
            failedSession = struct();
            failedSession.sessionName = sessionName;
            failedSession.message = extractionError.message;
            chamberMetadata.failedSessions(end + 1) = failedSession;
            summaryRows(end + 1, :) = {chamberName, sessionName, false, sessionOutputPath, extractionError.message}; %#ok<SAGROW>
        end
    end

    chamberOutputPath = fullfile(chamberOutputDir, ...
        sprintf('Behaviour_RAW_%s_AllSessions_%s.mat', chamberName, fileConditionTag));
    chamberMetadata.chamberOutputPath = chamberOutputPath;

    Data = struct();
    Data.Behav = chamberBehaviour;
    Data.meta = chamberMetadata;

    save(chamberOutputPath, 'Data', 'chamberBehaviour', 'chamberMetadata', '-v7.3');

    behaviourByChamber.(chamberName) = chamberBehaviour;
    metadataByChamber.(chamberName) = chamberMetadata;
end

extractionSummary = cell2table(summaryRows, ...
    'VariableNames', {'Chamber', 'Session', 'Success', 'OutputPath', 'ErrorMessage'});

allChambersOutputPath = fullfile(outputRoot, ...
    sprintf('Behaviour_RAW_AllChambers_%s.mat', fileConditionTag));

Data = struct();
Data.Behav = behaviourByChamber;
Data.meta = struct();
Data.meta.outputRoot = outputRoot;
Data.meta.chamberNames = chamberNames;
Data.meta.extractionPms = extractionPms;
Data.meta.conditionIndexMap = conditionIndexMap;
Data.meta.summary = extractionSummary;

save(allChambersOutputPath, ...
    'Data', 'behaviourByChamber', 'metadataByChamber', 'extractionSummary', '-v7.3');

fprintf('Behaviour extraction saved in %s\n', outputRoot);

function conditionIndexMap = buildBehaviourConditionIndexMap()
    behaviourIndex = (1:24)';
    conditionCode = [ones(8, 1); 2 * ones(8, 1); 3 * ones(8, 1)];
    targetDirection = repmat((1:8)', 3, 1);
    conditionName = strings(24, 1);
    activeMonkey = strings(24, 1);

    conditionName(1:8) = "Act S / Obs K";
    conditionName(9:16) = "Obs S / Act K";
    conditionName(17:24) = "Act S / Act K";

    activeMonkey(1:8) = "S";
    activeMonkey(9:16) = "K";
    activeMonkey(17:24) = "S and K";

    conditionIndexMap = table(behaviourIndex, conditionCode, targetDirection, ...
        conditionName, activeMonkey, ...
        'VariableNames', {'Index', 'ConditionCode', 'Direction', 'ConditionName', 'ActiveMonkey'});
end
