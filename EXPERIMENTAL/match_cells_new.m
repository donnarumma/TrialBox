addpath("GOODCELLS\")
%addpath("SpikesExtracted\")

load('GoodCells_K.mat')
load('GoodCells_S.mat')
%% PATH AN PATTERNS WHOM REFER TO
% next lines change according to path and data (name)
session_list=[22,25,33,35,36,38,42,43]
path_="C:\Users\loren\Desktop\USB_Content\AI_PhD_Neuro_CNR\Empirics\GIT_stuff\dati_Mirco\last_oct_25\SpikesExtracted"
cd(path_);
f_pattern = 'all_0_100_*.mat';

% MATCHING FILES
M_files = dir(f_pattern);
% patern of files to be loaded
%% extract and put in a structure the above specified sessions
all_sessions = cell(1, length(session_list));
for k = 1:length(M_files)
    fullname = fullfile(M_files(k).folder, M_files(k).name);

    fprintf('Loading %s\n', fullname);

    loaded = load(fullname);
    data_trials_ = loaded.data_trials;
    % Extract sessio number with a not found warning
    tokens = regexp(M_files(k).name, 'SK(\d+)H', 'tokens');
    if isempty(tokens)
        warning("File %s: can't find session number", M_files(k).name);
        continue;
    end
    sess_num = str2double(tokens{1}{1});

    idx = find(session_list == sess_num);

    if isempty(idx)
        warning('Session %d non è in session_list, salto.', sess_num);
        continue;
    end

  data_trials_ = loaded.data_trials;
  all_sessions{idx} = data_trials_;
end

%% find max and min values for condition and direction and select the 
% the appropriate neurons according to Goocells mat(s)
% output a new list of sessions including new fields (correct neurons for
% either monkey S and K)

all_cond=[]
all_dir=[]
for k=1:numel(session_list)
    all_cond= [all_cond, [all_sessions{k}.trialTypeCond]];
    all_dir = [all_dir, [all_sessions{k}.trialTypeDir]];
end
max_dir=max(all_dir)
min_dir=min(all_dir)
max_cond=max(all_cond)
min_cond=min(all_cond)
% tag direction-condition combo
for k=1:numel(session_list)
    curr_sess_id=session_list(k)
    curr_sess=all_sessions{k};
    all_n_id=curr_sess(1).All_N_Identity

    % matching neuroni
    Good_K=GoodCells_K(GoodCells_K(:,1)==curr_sess_id,:);
    Good_S=GoodCells_S(GoodCells_S(:,1)==curr_sess_id,:);
   
    [~,~,b_k]=intersect(Good_K(:,2:3),  all_n_id,'rows')
    [~,~,b_s]=intersect(Good_S(:,2:3),  all_n_id,'rows')

    N=numel(curr_sess)

    combo_counter = zeros(1, max_cond*max_dir); 
    for i = 1:N
        Spikes=curr_sess(i).Spikes;
        % 24 total combis
            
        cond=curr_sess(i).trialTypeCond
        dir=curr_sess(i).trialTypeDir
        tag=(cond-1)*max_dir+dir;
        combo_counter(tag)=combo_counter(tag)+1
        curr_sess(i).tag=combo_counter(tag);
        curr_sess(i).S_spikes=Spikes(b_s,:);
        curr_sess(i).K_spikes=Spikes(b_k,:);

    end
    all_sessions_new{k}=curr_sess
end
            

%% Generating S and K data matrix joining session
all_tags = [];
for j = 1:numel(all_sessions_new)
    all_tags = [all_tags, [all_sessions_new{j}.tag]];
end

max_tag = max(all_tags);
min_tag = min(all_tags);

S_final = cell(max_cond, max_dir, max_tag);
K_final = cell(max_cond, max_dir, max_tag);

for j = 1:numel(all_sessions_new)
    curr_sess = all_sessions_new{j};

    for i = 1:numel(curr_sess)
        cond = curr_sess(i).trialTypeCond;
        dir  = curr_sess(i).trialTypeDir;
        tag  = curr_sess(i).tag;
                S_final{cond, dir, tag} = [S_final{cond, dir, tag}; curr_sess(i).S_spikes];
                K_final{cond, dir, tag} = [K_final{cond, dir, tag}; curr_sess(i).K_spikes];
            end
        end
  
%%%%%% Matrices including Neural and Behavioral activity
S_matrix=[]
K_matrix=[]
   for t=min_tag:max_tag
       for c=min_cond:max_cond
           for d=min_dir:max_dir
               s_m=[]
               k_m=[]
               spikes_S=S_final{c,d,t}';
               ll_S=size(spikes_S,1);
               s_m=[spikes_S,ones(ll_S,1)*c ones(ll_S,1)*d];
               S_matrix=[S_matrix;s_m];
               

               spikes_K=K_final{c,d,t}';
               ll_K=size(spikes_K,1);
               k_m=[spikes_K,ones(ll_K,1)*c ones(ll_K,1)*d];
               K_matrix=[K_matrix;k_m];
             
           end
       end
   end
 
%% Change to suit personal preferences
save( 'dati_mirco_15_11_k_partial','K_matrix')
save( 'dati_mirco_15_11_s_partial','S_matrix')


%% S data arrangement to fit cebra
K=1
trial_id_s=ones(1,length(S_matrix))'
for i =2:length(S_matrix)

    if S_matrix(i-1,end)==S_matrix(i,end)
        trial_id_s(i)=K;
    else
        K=K+1
        trial_id_s(i)=K;
        
    end
end
%trial_id=trial_id'

S_matrix=[S_matrix trial_id_s];

S_matrix(:,end-1)
% assuming condition is declared in the third to last col
s_active_neural=S_matrix(S_matrix(:,end-2)==1,1:end-3)
% assume trial counter is the last column
s_active_trial_id=S_matrix(S_matrix(:,end-2)==1,end);
% assume trial value (i.e. label) is the 2nd to last col
s_active_trial=S_matrix(S_matrix(:,end-2)==1,end-1);

s_passive_neural=S_matrix(S_matrix(:,end-2)==2,1:end-3);
s_passive_trial_id=S_matrix(S_matrix(:,end-2)==2,end);
s_passive_trial=S_matrix(S_matrix(:,end-2)==2,end-1);

s_joint_neural=S_matrix(S_matrix(:,end-2)==3,1:end-3);
s_joint_trial_id=S_matrix(S_matrix(:,end-2)==3,1:end);
s_joint_trial=S_matrix(S_matrix(:,end-2)==3,end-1);

% Change to suit personal preferences
save('dati_mirco_15_11_400_s',"s_active_neural","s_active_trial_id","s_active_trial", ...
    "s_passive_neural","s_passive_trial_id","s_passive_trial", ...
    "s_joint_neural","s_joint_trial_id","s_joint_trial")

%% K data arrangement to fit cebra
K=1
trial_id_k=ones(1,length(K_matrix))'
for i =2:length(K_matrix)

    if K_matrix(i-1,end)==K_matrix(i,end)
        trial_id_k(i)=K;
    else
        K=K+1
        trial_id_k(i)=K;
        
    end
end
%trial_id=trial_id'

K_matrix=[K_matrix trial_id_k];

k_passive_neural=K_matrix(K_matrix(:,end-2)==1,1:end-3)
k_passive_trial_id=K_matrix(K_matrix(:,end-2)==1,end);
k_passive_trial=K_matrix(K_matrix(:,end-2)==1,end-1);

k_active_neural=K_matrix(K_matrix(:,end-2)==2,1:end-3);
k_active_trial_id=K_matrix(K_matrix(:,end-2)==2,end);
k_active_trial=K_matrix(K_matrix(:,end-2)==2,end-1);

k_joint_neural=K_matrix(K_matrix(:,end-2)==3,1:end-3);
k_joint_trial_id=K_matrix(K_matrix(:,end-2)==3,1:end);
k_joint_trial=K_matrix(K_matrix(:,end-2)==3,end-1);

% Change to suit personal preferences

save('dati_mirco_15_11_400_k',"k_active_neural","k_active_trial_id","k_active_trial", ...
    "k_passive_neural","k_passive_trial_id","k_passive_trial", ...
    "k_joint_neural","k_joint_trial_id","k_joint_trial")


