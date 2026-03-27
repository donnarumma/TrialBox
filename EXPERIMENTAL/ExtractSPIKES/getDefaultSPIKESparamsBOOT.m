function params = getDefaultSPIKESparamsBOOT()
% function params = getDefaultSPIKESparamsBOOT()

VARIABLES_monkey_joint_action;

%% --- DEFAULT SELECTION ---
params.chamber = 'Frontal';   % 'Frontal' | 'Parietal'

params.interval        =[-200,500]; % time interval 
                                    % Time max overall 3.1 s or 2.2

params.selMonkey = 'all'; % selection Monkey "S" or "K"
                          % "all" both S and K
%% --- SESSION LISTS ---
sessions.Frontal = { ...
                    'SK007','SK008','SK009','SK010','SK011','SK012','SK013','SK014', ...
                    'SK020','SK021','SK022','SK023','SK024','SK025','SK027','SK028', ...
                    'SK029','SK030','SK031','SK032','SK033','SK034','SK035','SK036', ...
                    'SK037','SK038','SK039','SK040','SK041','SK042','SK043','SK044', ...
                    'SK045','SK046','SK047','SK048' };

sessions.Parietal = { ...
                    'SK003','SK004','SK049','SK050','SK051','SK052','SK053','SK054', ...
                    'SK055','SK056','SK057','SK058','SK059','SK060','SK061','SK062', ...
                    'SK063','SK064','SK065','SK066','SK067','SK069','SK070','SK071', ...
                    'SK072' };

%% --- SELECT SESSION LIST BASED ON CHAMBER ---
if ~isfield(sessions, params.chamber)
    error('Unknown Chamber "%s". Allowed: %s', ...
        params.chamber, strjoin(fieldnames(sessions), ', '));
end

params.session_list = sessions;
%% --- ALIGNMENT MODE ---
params.MODE_ALIGN = TARGET_APPEARS_MODE; 
% or MOVEMENT_ONSET_MODE; 
% MOVEMENT_ONSET_K_MODE;
% MOVEMENT_ONSET_S_MODE
%% SAVE PARAMS
params.singleDir_save = 'SingleSession';
params.dir_save = 'AllSessions';
