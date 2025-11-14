addpath("GOODCELLS\")
addpath("SpikesExtracted\")
%% Questa parte va automatizata (anonimizzata)
load('GoodCells_K.mat')
load('GoodCells_S.mat')
all_trials_22=load('SpikesExtracted\all_200_1000_frontalSK022H_SUA.mat')
all_trials_33=load('SpikesExtracted\all_200_1000_frontalSK033H_SUA.mat')
all_trials_25=load('SpikesExtracted\all_200_1000_frontalSK025H_SUA.mat')
all_trials_35=load('SpikesExtracted\all_200_1000_frontalSK035H_SUA.mat')
all_trials_36=load('SpikesExtracted\all_200_1000_frontalSK036H_SUA.mat')
all_trials_38=load('SpikesExtracted\all_200_1000_frontalSK038H_SUA.mat')
all_trials_42=load('SpikesExtracted\all_200_1000_frontalSK042H_SUA.mat')
all_trials_43=load('SpikesExtracted\all_200_1000_frontalSK043H_SUA.mat')

all_session={all_trials_22.data_trials,all_trials_25.data_trials, ...
    all_trials_33.data_trials,all_trials_35.data_trials, ...
    all_trials_36.data_trials,all_trials_38.data_trials, ...
    all_trials_42.data_trials,all_trials_43.data_trials
    }

session_list=[22,25,33,35,36,38,42,43]
%%
all_session_new={}
all_dir=[]
all_cond=[]

for k=1:numel(session_list)
    all_cond= [all_cond, [all_session{k}.trialTypeCond]];
    all_dir = [all_dir, [all_session{k}.trialTypeDir]];
end
max_dir=max(all_dir)
min_dir=min(all_dir)
max_cond=max(all_cond)
min_cond=min(all_cond)
% tagghiamo le combo direzione condizione
for k=1:numel(session_list)
    curr_sess_id=session_list(k)
    curr_sess=all_session{k};
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
        % 24 combo totali
            
        cond=curr_sess(i).trialTypeCond
        dir=curr_sess(i).trialTypeDir
        tag=(cond-1)*max_dir+dir;
        combo_counter(tag)=combo_counter(tag)+1
        curr_sess(i).tag=combo_counter(tag);
        curr_sess(i).S_spikes=Spikes(b_s,:);
        curr_sess(i).K_spikes=Spikes(b_k,:);

    end
    all_session_new{k}=curr_sess
end
            

%% Genero dati S e K unendo le sessioni new
% come prima per cond e dir contare il max e il min tag

all_tags = [];
for j = 1:numel(all_session_new)
    all_tags = [all_tags, [all_session_new{j}.tag]];
end

max_tag = max(all_tags);
min_tag = min(all_tags);

S_final = cell(max_cond, max_dir, max_tag);
K_final = cell(max_cond, max_dir, max_tag);

for j = 1:numel(all_session_new)
    curr_sess = all_session_new{j};

    for i = 1:numel(curr_sess)
        cond = curr_sess(i).trialTypeCond;
        dir  = curr_sess(i).trialTypeDir;
        tag  = curr_sess(i).tag;
                S_final{cond, dir, tag} = [S_final{cond, dir, tag}; curr_sess(i).S_spikes];
                K_final{cond, dir, tag} = [K_final{cond, dir, tag}; curr_sess(i).K_spikes];
            end
        end
  
%%%%%% matrici finali Attività neurale+Comportamento
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
 




