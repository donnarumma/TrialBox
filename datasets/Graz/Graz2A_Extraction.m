function   Graz2A_Extraction(infopath)
% function Graz2A_Extraction(infopath)
% Dataset eeegplanetion: "Leeb, R., Brunner, C., Müller-Putz, G., Schlögl, A., & Pfurtscheller, G. J. G. U. O. T. (2008). BCI Competition 2008–Graz data set B.
%       Graz University of Technology, Austria, 16, 1-6."

% Link for data:    www.bbci.de/competition/iv/ -> Download of data sets -> agree submit -> Data sets 2a: ‹4-class motor imagery>
sourcepath=infopath.source; % where are the gdf files (training and test data)

% Link for labels:  www.bbci.de/competition/iv/ -> News 
%                   -> Results of the BCI Competition IV 

labelpath =infopath.labels;  % label mat files

%                   -> True Labels of Competition's Evaluation Sets -> Data sets 2a: 
% Link to BioSig:   -> https://biosig.sourceforge.net/download.html
% addpath -> path to ('biosig/t200_FileAccess/');
% addpath -> path to ('biosig/t250_ArtifactPreProcessingQualityControl/');

savepath  =infopath.save;    % output mat files 

% Specifics:
% 1:    Left Hand
% 2:    Right Hand
% 3:    Foot
% 4:    Tongue

%% Extraction data
% define path of the data extracted by link
% sourcepath  = 'D:\D_Ausilio EEG\Confronto Graz\dataGraz\';
% sourcepath  = '~/DATA/GRAZ/BCICIV_2a_gdf/';

% define labels true path
% labelpath   = 'D:\D_Ausilio EEG\Confronto Graz\true_labels\';
% savepath    = 'D:\D_Ausilio EEG\Confronto Graz\Extracted_Data\';
% savepath    = '~/DATA/GRAZ/BCICIV_2a_mat/';
% labelpath   = '~/DATA/GRAZ/true_labels/';
% load training file T and test file E (Evaluation)
file = dir(sourcepath);

T_pattern = '^A[0][1-9]T.*\.gdf$';
E_pattern = '^A[0][1-9]E.*\.gdf$';

file_T = struct();
file_E = struct();
nT = 1;
nE = 1;
for i = 1:length(file)
    if ~isempty(regexp(file(i).name, T_pattern, 'once'))
       file_T(nT).name = file(i).name;
       nT = nT+1;
    elseif ~isempty(regexp(file(i).name, E_pattern, 'once'))
       file_E(nE).name = file(i).name;
       nE = nE+1;
    end
end

% time interval with: 0 in Visual Cue, tstart at FixationCross Tend at the end of
% Motor Imagery task
tstart  = -2;
tend    =  6;

for idsubj = 1:length(file_T)
    cl=load(fullfile(labelpath,replace(file_T(idsubj).name,'.gdf','.mat')));
    classlabel=cl.classlabel;
    % Train data
    fprintf('Extracting %s\n',file_T(idsubj).name);
    % [sT, hT]    = sload(file_T(idsubj).name);
    [sT, hT]    = sload([sourcepath file_T(idsubj).name]);
    fc          = hT.SampleRate;
    nTr = 1;
    EEG_train = struct();
    for nev_tr = 1:length(hT.EVENT.TYP)
        timeclass_tr = hT.EVENT.POS(nev_tr);
        if hT.EVENT.TYP(nev_tr)==769 % Class 1 Left Hand
            EEG_train(nTr).eeg = Graz2A_getData(sT,timeclass_tr,tstart,tend,fc);
            EEG_train(nTr).label = 1;
            nTr = nTr +1;
        elseif hT.EVENT.TYP(nev_tr)==770 % Class 2 Rigth Hand
            EEG_train(nTr).eeg = Graz2A_getData(sT,timeclass_tr,tstart,tend,fc);
            EEG_train(nTr).label = 2;
            nTr = nTr +1;
        elseif hT.EVENT.TYP(nev_tr)==771 % Class 3 Foot
            EEG_train(nTr).eeg = Graz2A_getData(sT,timeclass_tr,tstart,tend,fc);
            EEG_train(nTr).label = 3;
            nTr = nTr +1;
        elseif hT.EVENT.TYP(nev_tr)==772 % Class 4 Tongue
            EEG_train(nTr).eeg = Graz2A_getData(sT,timeclass_tr,tstart,tend,fc);
            EEG_train(nTr).label = 4;
            nTr = nTr +1;
        end
    end

    for iTr = 1:length(EEG_train)
        EEG_train(iTr).eeg = cell2mat(EEG_train(iTr).eeg);
        if EEG_train(iTr).label ~= classlabel(iTr)
            disp('error')
            pause
        end
    end

    % Test Data
    fprintf('Extracting %s...',file_E(idsubj).name);
    [sE, hE]            = sload([sourcepath file_E(idsubj).name]);
    fc                  = hE.SampleRate;
    nTs = 1;
    EEG_test = struct();
    for nev_ts = 1:length(hE.EVENT.TYP)
        timeclass_ts = hE.EVENT.POS(nev_ts);
        if hE.EVENT.TYP(nev_ts)==783 % Unknown class only Visual Cue event
            EEG_test(nTs).eeg = Graz2A_getData(sE,timeclass_ts,tstart,tend,fc);
            nTs = nTs +1;
        end
    end
    for iTr = 1:length(EEG_test)
        EEG_test(iTr).eeg = cell2mat(EEG_test(iTr).eeg);
        EEG_test(iTr).truelab = classlabel(iTr);
    end
    data_raw.EEG_train = EEG_train;
    data_raw.EEG_test = EEG_test;
    data_raw.fsample = fc;

    % save the data for the further analysis
    filename = replace(file_T(idsubj).name,'.gdf','E.mat');
    fullpath = fullfile(savepath, filename);
    fprintf('Saving %s...\n',fullpath);
    save(fullpath, 'data_raw');
end
