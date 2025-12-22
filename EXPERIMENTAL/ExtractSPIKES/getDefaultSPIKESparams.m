function params = getDefaultSPIKESparams()
% function params = getDefaultSPIKESparams()

VARIABLES_monkey_joint_action;
params.interval        =[-199,400]; % time interval 
                                    % Time max overall 3.1 s or 2.2

params.selMonkey = 'all'; % selection Monkey "S" or "K"
                          % "all" both S and K
% All frontal Session
params.session_list = {'SK011','SK012','SK021','SK022','SK023',...
                        'SK024','SK025','SK027','SK028','SK030','SK032',...
                        'SK033','SK034','SK035','SK036','SK037','SK038',...
                        'SK039','SK040','SK041','SK042','SK043','SK044',...
                        'SK045','SK046','SK047','SK048'};

% Frontal Sessions with only 192 Trials 
% params.session_list = {'SK022','SK033','SK025','SK035','SK036',...
%                         'SK038','SK042','SK043'};

params.MODE_ALIGN = TARGET_APPEARS_MODE; 
% or MOVEMENT_ONSET_MODE; MOVEMENT_ONSET_MODE; MOVEMENT_ONSET_K_MODE;
% MOVEMENT_ONSET_S_MODE
