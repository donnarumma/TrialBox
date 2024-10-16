function   par = BATTAGLIA_CSDParams(parSource)
% function par = BATTAGLIA_CSDParams(parSource)

par.ifserver    = true;     % do not plot result on server
par.whichmodel  = 7;        % 7 Model for action session. 5 Model for all sessions
par.isdemo      = 1;        % getSelectionIndexes(1,1:3);
par.dmode       = 7;        % LFP extraction mode 7
par.selS        = 1;        % LFP index choice
par.selK        = [];
par.session_name= 'SK009';  % which session
par.custom_model= [];
par.donlfp      = false;

try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
