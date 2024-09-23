function   [pE,pC] = customPriorsv3(pE,pC,par)
% function [pE,pC] = customPriorsv3(pE,pC,TEST_DIR)
% try
%     TEST_DIR;
% catch
%     TEST_DIR = '~/TESTS/SAPIENZA/DCM/';
% end
%S_PRIORS='LFP_M7_SK009_D7_S1_K0_C_1';
%S_PRIORS    =par.S_PRIORS;
%DC_load     =load([TEST_DIR S_PRIORS]);
DC_load     =load(par.S_PRIORS);
DCM_S       =DC_load.DCM;

% K_PRIORS='LFP_M7_SK009_D7_S0_K1_C_9';
% K_PRIORS    =par.K_PRIORS;
% DC_load     =load([TEST_DIR K_PRIORS]);
DC_load     =load(par.K_PRIORS);
DCM_K       =DC_load.DCM;

EpS=DCM_S.Ep;
CpS=DCM_S.Pp;
EpK=DCM_K.Ep;
CpK=DCM_K.Pp;

Ep.R = mean([EpS.R;EpK.R]);
Ep.T = [EpS.T;EpK.T];
Ep.G = [EpS.G;EpK.G];
Ep.H = [EpS.H;EpK.H];
Ep.A = pE.A;
for iA=1:length(pE.A)
    Ep.A{iA}(1,1)=EpS.A{iA};
    Ep.A{iA}(2,2)=EpK.A{iA};
end
%% WARNING
Ep.B = pE.B;
idir = 1;
for iB=1:length(pE.B)
    Ep.B{iB}(1,1) = EpS.B{idir};
    Ep.B{iB}(2,2) = EpK.B{idir};
end
fprintf('WARNING\n');
Ep.C    = pE.C;
Ep.D    = pE.D;
Ep.I    = [EpS.I + EpK.I]/2;
Ep.L    = [EpS.L,EpK.L];
Ep.a    = [EpS.a , EpK.a];
Ep.b    = [EpS.b + EpK.b]/2;
Ep.c    = [EpS.c + EpK.c]/2;
Ep.d    = [EpS.d , EpK.d];
Ep.f    = [EpS.f + EpK.f]/2;
Ep.Lpos = [EpS.Lpos + EpK.Lpos]/2;
Ep.J    = [EpS.J + EpK.J]/2;

%
pE   =Ep;
pC   =pC;












