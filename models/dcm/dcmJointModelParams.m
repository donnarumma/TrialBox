function   par = dcmJointModelParams(parSource)
% function par = BattagliaArrangeTrialsParams(parSource)

par.mstep       = 8;        % M-step: Fisher scoring scheme

try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end

