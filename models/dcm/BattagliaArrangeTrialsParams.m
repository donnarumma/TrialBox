function   par = BattagliaArrangeTrialsParams(parSource)
% function par = BattagliaArrangeTrialsParams(parSource)

% par.whichmodel  = 7;        % 7 Model for action session. 5 Model for all sessions
par.isdemo      = 1;        % getSelectionIndexes(1,1:3);
par.dmode       = 7;        % LFP extraction mode 7
par.selS        = 1;        % LFP index choice monkey S
par.selK        = 1;        % LFP index choice monkey K;
par.session_name= [];       % which session (e.g.'SK009')

try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end

