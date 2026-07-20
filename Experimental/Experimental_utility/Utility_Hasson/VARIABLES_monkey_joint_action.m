%function VARIABLES_monkey_joint_action
SOLO_S_CONDITION=8001;
SOLO_K_CONDITION=8002;
T1              =8101;
T2              =8102;
T3              =8103;
T4              =8104;
T5              =8105;
T6              =8106;
T7              =8107;
T8              =8108;
T_LABELS        ={'T1','T2','T3','T4','T5','T6','T7','T8'};
C_LABELS        ={'JOINT','SOLO S','SOLO K'};
JOINT_CONDITION =8004;
SUCCESS         =3000;
TARGET_APPEARS  =2024;
S_MOVEMENT_ONSET=2107;
K_MOVEMENT_ONSET=2207;
CHAMBERS        ={'frontal','parietal'};
ET1S 	 	=2206; % entrata nel target periferico da parte del cursore 1 (scimmia S)
ET2K		=2106; % entrata nel target periferico da parte del cursore 2 (scimmia K)

% forse il leader scemo è S - no è K
% analizzare come le componenti principali di S spiegano K e viceversa.
% ricostruire S o K a partire dalle componenti principali S+K
MOVEMENT_ONSET_MODE            = 0;
TARGET_APPEARS_MODE            = 1;
MOVEMENT_ONSET_K_MODE          = 2;
MOVEMENT_ONSET_S_MODE          = 3;
GOOD_CELL_SELECTION            = 4;
ALL_CELL_SELECTION             = 5;
TUNING_CELL_SELECTION          = 6;
TUNING_CELL_SOLO_SELECTION     = 7;
TUNING_CELL_OBS_SELECTION      = 8;
TUNING_CELL_JOINT_SELECTION    = 9;
TUNING_CELL_NO_SOLO_SELECTION  =10;
TUNING_CELL_NO_OBS_SELECTION   =11;
TUNING_CELL_NO_JOINT_SELECTION =12;
TUNING_CELL_NO_SELECTION       =13;
TUNING_RT                      =14;
TUNING_MT                      =15;
TUNING_RT_MT                   =16;
