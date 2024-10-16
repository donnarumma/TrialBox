function   out=csdSessionExtract(session_name,dmode,selS,selK)
% function out=csdSessionExtract(session_name,dmode,selS,selK)
% example:
% session_name = 'SK009';
try
    dmode;
catch
    dmode = 7;
end
params       = getLFPparamsDonnarumma(dmode);
behav_JOINT  = 'H';
behav_EYE    = 'E';
signal_type  = 'Raw';
f_Joint      = load([session_name behav_JOINT '_' signal_type]);
Raw_H        = f_Joint.Trials;

try 
    f_Eye        = load([session_name behav_EYE   '_' signal_type]);
    Raw_E        = f_Eye.Trials;
catch
    Raw_E        = [];
end
Raw_H        = EXTRACT_LFP(Raw_H, Raw_E, params);

tstart       = -0.7;
tend         =  0.7;
TA           = 2024; % align at Target Appears (periferic)
Trials       = Donnarumma_Align_Trials(Raw_H, TA, tstart,tend);

fprintf('Extracting session %s - Chamber: %s\n',session_name,Trials(1).Chamber);

%% select the good LFP indexes
Depth        = Trials(1).Depth;
iLFP         = false(1,length(Depth));
dLFP         = ~isnan(Depth);

vLFP_S       = iLFP;
vLFP_S(1:5)  = true;
vLFP_S       = vLFP_S & dLFP;
iLFP_S       = find(vLFP_S);

%fprintf('iLFP S\n'); disp(iLFP_S);

vLFP_K       = iLFP;
vLFP_K(6:10) = true;
vLFP_K       = vLFP_K & dLFP;
iLFP_K       = find(vLFP_K);

%fprintf('iLFP K\n'); disp(iLFP_K);
try
    selIndex     = [iLFP_S(selS),iLFP_K(selK)];   % take the first good LFP for each monkey
catch
    out          = [];
    return
end
%% LABELS

[condition, direction, klabels, success] = getClassInfo(Trials);
Trials       = Trials(success);
condition    = condition(success,:);
direction    = direction(success,:);
klabels      = klabels(success,:);

NTrials      = length(Trials);
dt           = mean(diff(Trials(1).RawTime)); % sample time

ARp          = 8; 
rsfactor     = 0.2;
Hz           = 1:64;

% CSD          = cell(1,NTrials);
LFP          = cell(1,NTrials);
t=tic;
for indTrial     = 1:NTrials
    % fprintf('Selecting Trial %g\n',indTrial);
    LFP{indTrial} = Trials(indTrial).LFP(selIndex,:)'; % time length x number of sources
    % CSD{indTrial} = DONNARUMMA_estimateCSD(LFP{indTrial},Hz,dt,ARp,rsfactor);
end
CSD=computeCSD_fromLFP(LFP,Hz,dt,ARp,rsfactor);

fprintf('Computed Cross Spectral Densities in %g s\n',toc(t));
out.CSD         = CSD;
out.LFP         = LFP;
out.Trials      = Trials;
out.condition   = condition;
out.direction   = direction;
out.klabels     = klabels;
out.dt          = dt;
out.Hz          = Hz;
out.selIndex    = selIndex;
out.ARp         = ARp;
out.rsfactor    = rsfactor;
out.iLFP_K      = iLFP_K;
out.iLFP_S      = iLFP_S;