%% MAIN_ExtractAllSPIKES
% copia della funzione S_translation di Donnarumma
% LISTA DEI PATH NECESSARI:
% PATH contenente i dati DATI SUA
% CARTELLA ExtractSPIKES in TrialBox

clear;
close all

[~,nn]=system('hostname'); nn=strtrim(nn);
if strcmp(nn,'rk018610')
    isonserver=true;
else
    isonserver=false;
end

params      = getDefaultSPIKESparams;

interval        = params.interval; 
selMonkey       = params.selMonkey;
session_list    = params.session_list;
MODE_ALIGN      = params.MODE_ALIGN;

% Name of Saved data
name        = 'AllFrontal';


VARIABLES_monkey_joint_action;
for n_sess = 1:length(session_list)
    session_name = session_list{n_sess};
    try
        session_name = strcat(session_name,'H_SUA');
    catch
        fprintf('Spikes_translation Warning: using %s\n',session_name);
        pause;
    end
    d       = load(session_name);
    Trials  = d.Trials;
    chamber = Trials(1).Chamber;
    tic;
    fprintf('Getting data from %s\n',session_name);

    DataSpike = getSpikesJointAllNeurons(Trials);
    fprintf('Elapsed time: %g s\n',toc(tic))
    DataSpike = getSuccessfullTrials(DataSpike,Trials);
    DataSpike = getJointConditionTrials(DataSpike,Trials);
    DataSpike = getSolo_S_ConditionTrials(DataSpike,Trials);
    DataSpike = getSolo_K_ConditionTrials(DataSpike,Trials);
    DataSpike = getAreaTrials(DataSpike,Trials);
    DataSpike = getGoodNeurons(DataSpike,Trials);
    % DataSpike=      getTuningSoloNeurons(S,Trials);
    % DataSpike=     getTuningJointNeurons(S,Trials);
    % DataSpike=       getTuningObsNeurons(S,Trials);
    % DataSpike=    getTuningNoSoloNeurons(S,Trials);
    % DataSpike=   getTuningNoJointNeurons(S,Trials);
    % DataSpike=     getTuningNoObsNeurons(S,Trials);
    % if length(CELL_SELECTION)>1
    %     DataSpike=          getTuningNeurons(S,Trials,CELL_SELECTION(2));
    %     DataSpike=        getNoTuningNeurons(S,Trials,CELL_SELECTION(2));
    % end
    DataSpike = getTimeMovementOnsetTargetTrials(DataSpike,Trials); %% new entry
    for ind=1:8
        DataSpike = getTConditionTrials(DataSpike,Trials,ind);
    end
    DataSpike = getTimeTargetTrials(DataSpike,Trials);
    %S=retuneSpikeWithTimeStamp1(S,Trials);
    if MODE_ALIGN==TARGET_APPEARS_MODE
        DataSpike = retuneSpikeWithTargetAppears(DataSpike,Trials);
    elseif MODE_ALIGN==MOVEMENT_ONSET_MODE
        DataSpike = retuneSpikeWithMovementOnset(DataSpike,Trials);
    elseif MODE_ALIGN==MOVEMENT_ONSET_K_MODE
        DataSpike = retuneSpikeWithMovementOnsetK(DataSpike,Trials);
    elseif MODE_ALIGN==MOVEMENT_ONSET_S_MODE
        DataSpike = retuneSpikeWithMovementOnsetS(DataSpike,Trials);
    end

    DataSpike = prepareSpikeAndTrials(DataSpike,Trials,interval,selMonkey);

    dat = DataSpike.dat;
    data_session = dat;

    % time vector
    t = linspace(interval(1), interval(2), size(dat(1).spikes,2));

    dataSession_trials = struct();

    for iTr=1:length(data_session)
        dataSession_trials(iTr).trialId                = data_session(iTr).trialId;
        dataSession_trials(iTr).timeSpikes             = t;
        dataSession_trials(iTr).Spikes                 = data_session(iTr).spikes;
        dataSession_trials(iTr).trialTypeDir           = data_session(iTr).trialType;
        dataSession_trials(iTr).trialTypeCond          = data_session(iTr).trialType2;
        dataSession_trials(iTr).SessName               = session_list{n_sess};
        dataSession_trials(iTr).chamber                = chamber;
    end

    dataSession_trials(1).All_N_Identity               = DataSpike.AllNeuronIdentity;
    dataSession_trials(1).All_Good_N_Identity          = DataSpike.AllGoodNeuronIdentity;

    dataSession_trials                                 = conditionalligne(dataSession_trials);
    % saving single session data
    save_path = fullfile(pwd, [name,filesep,'SpikesExtracted']);

    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end

    name_save = strcat(selMonkey,'_',num2str(interval(1)),'_', ...
        num2str(interval(2)),'_',DataSpike.area, session_name);
    save(fullfile(save_path, strcat(name_save,'.mat')), 'dataSession_trials');
end

f_pattern = [save_path, filesep,'all_',num2str(interval(1)),'_',num2str(interval(2)),'_frontal*.mat'];

% Full path
M_files = dir(f_pattern);

disp("Files found:");
disp({M_files.name});

% MATCHING FILES
M_files = dir(f_pattern);

sess_temp = erase(session_list, 'SK');
session_num = str2double(sess_temp);

dataessions = cell(1, length(session_num));
for k = 1:length(M_files)
    fullname = fullfile(M_files(k).folder, M_files(k).name);

    fprintf('Loading %s\n', fullname);

    loaded = load(fullname);
    data_trials_load = loaded.dataSession_trials;
    % Extract session number with a not found warning
    tokens = regexp(M_files(k).name, 'SK(\d+)H', 'tokens');
    if isempty(tokens)
        warning("File %s: can't find session number", M_files(k).name);
        continue;
    end
    sess_nb = str2double(tokens{1}{1});

    idx = find(session_num == sess_nb);

    if isempty(idx)
        warning('Session %d is not in session_list, skipping.', sess_nb);
        continue;
    end

    data_trials_load = loaded.dataSession_trials;
    dataessions{idx} = data_trials_load;
end

%% find max and min values for condition and direction and select the
% the appropriate neurons according to Goocells mat(s)
% output a new list of sessions including new fields (correct neurons for
% either monkey S and K)
totTrials = 0;
for k = 1:numel(session_num)
    totTrials = totTrials + numel(dataessions{k});
end

all_cond = zeros(1, totTrials);
all_dir  = zeros(1, totTrials);

idx = 1;

for k = 1:numel(session_num)

    curr_sess = dataessions{k};

    if isempty(curr_sess)
        error("Empty session %d: dataessions{%d} = []", session_num(k), k);
    end

    all_n_id = curr_sess(1).All_N_Identity;
    if size(unique(all_n_id,'rows'),1) ~= size(all_n_id,1)
        error("Duplicates in neurons! Sessione %d", session_num(k));
    end

    n = numel(curr_sess);

    all_cond(idx:idx+n-1) = [curr_sess.trialTypeCond];
    all_dir(idx:idx+n-1)  = [curr_sess.trialTypeDir];

    idx = idx + n;
end

max_dir     = max(all_dir);
min_dir     = min(all_dir);
max_cond    = max(all_cond);
min_cond    = min(all_cond);
% tag direction-condition
dataessions_new = cell(size(dataessions));
load('GoodCells_K.mat')
load('GoodCells_S.mat')
for k=1:numel(session_num)
    curr_sess_id    = session_num(k);
    curr_sess       = dataessions{k};
    all_n_id        = curr_sess(1).All_N_Identity;

    % matching neurons
    Good_K          = GoodCells_K(GoodCells_K(:,1)==curr_sess_id,:);
    Good_S          = GoodCells_S(GoodCells_S(:,1)==curr_sess_id,:);

    [~,~,intGood_k]       = intersect(Good_K(:,2:3),  all_n_id,'rows');
    [~,~,intGood_s]       = intersect(Good_S(:,2:3),  all_n_id,'rows');

    N               = numel(curr_sess);
    combo_counter   = zeros(1, max_cond*max_dir);

    for i = 1:N
        Spikes                                      = curr_sess(i).Spikes;
        cond                                        = curr_sess(i).trialTypeCond;
        dir                                         = curr_sess(i).trialTypeDir;
        tag                                         = (cond-1)*max_dir+dir;
        combo_counter(tag)                          = combo_counter(tag)+1;
        curr_sess(i).tag                            = combo_counter(tag);
        curr_sess(i).S_spikes                       = Spikes(intGood_s,:);
        curr_sess(i).K_spikes                       = Spikes(intGood_k,:);
        curr_sess(i).Neuron_S.All_N_Id_Position     = intGood_s;
        curr_sess(i).Neuron_S.Good                  = Good_S;
        curr_sess(i).Neuron_K.All_N_Id_Position     = intGood_k;
        curr_sess(i).Neuron_K.Good                  = Good_K;
    end

    dataessions_new{k}  = curr_sess;
end


%% Generating S and K data matrix joining session
totN = 0;
for j = 1:numel(dataessions_new)
    totN = totN + numel(dataessions_new{j});
end

all_tags = zeros(1, totN);   % oppure false / categorical

idx = 1;
for j = 1:numel(dataessions_new)
    curr                    = dataessions_new{j};
    n                       = numel(curr);
    all_tags(idx:idx+n-1)   = [curr.tag];
    idx                     = idx + n;
end

max_tag = max(all_tags);
min_tag = min(all_tags);

S_final = cell(max_cond, max_dir, max_tag);
K_final = cell(max_cond, max_dir, max_tag);

trialId_final           = cell(max_cond, max_dir, max_tag);
trialSession_final      = cell(max_cond, max_dir, max_tag);
timeSpikes_final        = cell(max_cond, max_dir, max_tag);
chamber_final           = cell(max_cond, max_dir, max_tag);
Neuron_final_S          = cell(max_cond, max_dir, max_tag);
Neuron_final_K          = cell(max_cond, max_dir, max_tag);

for j = 1:numel(dataessions_new)
    curr_sess = dataessions_new{j};
    for i = 1:numel(curr_sess)
        cond = curr_sess(i).trialTypeCond;
        dir  = curr_sess(i).trialTypeDir;
        tag  = curr_sess(i).tag;
        
        S_final{cond, dir, tag}             = [S_final{cond, dir, tag}; curr_sess(i).S_spikes];
        K_final{cond, dir, tag}             = [K_final{cond, dir, tag}; curr_sess(i).K_spikes];
        trialId_final{cond, dir, tag}       = [trialId_final{cond, dir, tag}; curr_sess(i).trialId];
        trialSession_final{cond, dir, tag}  = [trialSession_final{cond, dir, tag}; {curr_sess(i).SessName}];
        timeSpikes_final{cond, dir, tag}    = [timeSpikes_final{cond, dir, tag}; curr_sess(i).timeSpikes];
        chamber_final{cond, dir, tag}       = [chamber_final{cond, dir, tag};{curr_sess(i).chamber}];
        Neuron_final_S{cond, dir, tag}      = [Neuron_final_S{cond, dir, tag};curr_sess(i).Neuron_S];
        Neuron_final_K{cond, dir, tag}      = [Neuron_final_K{cond, dir, tag};curr_sess(i).Neuron_K];
    end
end

data = struct( ...
    'trialID',         [], ... 
    'trialIdSession',  {}, ...  
    'timeSpikes',      [], ...
    'trialTypeCond',   [], ...
    'trialTypeDir',    [], ...
    'S_Spikes',        [], ...
    'K_Spikes',        [], ...
    'Tag',             []);

prod = 0;

for cc = 1:max_cond
    for dd = 1:max_dir
        for tt = 1:max_tag
            S_s = S_final{cc,dd,tt};
            K_k = K_final{cc,dd,tt};

            if isempty(S_s) && isempty(K_k)
                continue
            end

            prod = prod + 1;

            data(prod).trialID          = trialId_final{cc,dd,tt};
            data(prod).trialIdSession   = trialSession_final{cc,dd,tt};
            data(prod).timeSpikes       = timeSpikes_final{cc,dd,tt};
            data(prod).trialTypeCond    = cc;
            data(prod).trialTypeDir     = dd;
            data(prod).Tag              = tt;

            data(prod).S_Spikes         = S_s;
            data(prod).K_Spikes         = K_k;

            data(prod).chamber          = chamber_final{cc,dd,tt};
            data(prod).Neuron.S         = Neuron_final_S{cc,dd,tt};
            data(prod).Neuron.K         = Neuron_final_K{cc,dd,tt};

        end
    end
end

%% tag filtering
data_trials = data([data.Tag]<=8);

save([name,'.mat'],'data_trials')