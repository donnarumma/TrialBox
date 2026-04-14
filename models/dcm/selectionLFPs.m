function [selS,selK,depthS,depthK] = selectionLFPs(session_name,par)

behav_JOINT = 'H';
signal_type = 'Raw';
method = lower(par.method);

f_Joint = load([session_name behav_JOINT '_' signal_type]);
Raw_H   = f_Joint.Trials;

%% Select good LFP indexes
Depth = Raw_H(1).Depth;
dLFP  = ~isnan(Depth);

vLFP_S = false(1,length(Depth));
vLFP_S(1:5) = true;
vLFP_S = vLFP_S & dLFP;

vLFP_K = false(1,length(Depth));
vLFP_K(6:10) = true;
vLFP_K = vLFP_K & dLFP;

validIdx_S   = find(vLFP_S);
validIdx_K   = find(vLFP_K);
validDepth_S = Depth(validIdx_S);
validDepth_K = Depth(validIdx_K);

if strcmp(method,'maximum')

    [depthS, selS] = max(validDepth_S);
    [depthK, selK] = max(validDepth_K);

elseif strcmpi(par.method,'minimum')

    [depthS, selS] = min(validDepth_S);
    [depthK, selK] = min(validDepth_K);

elseif strcmp(method,'mean')

    [selS, depthS] = selectClosestToTarget(validDepth_S, mean(validDepth_S));
    [selK, depthK] = selectClosestToTarget(validDepth_K, mean(validDepth_K));

elseif strcmp(method,'median')

    [selS, depthS] = selectClosestToTarget(validDepth_S, median(validDepth_S));
    [selK, depthK] = selectClosestToTarget(validDepth_K, median(validDepth_K));

elseif strcmp(method,'all')

    depthS = validDepth_S;
    depthK = validDepth_K;
    selS   = 1:length(validIdx_S);
    selK   = 1:length(validIdx_K);

else
    error('selectionLFPs:UnknownMethod', ...
        'Metodo non riconosciuto: %s', par.method);
end

end

function [selRel, depthSel] = selectClosestToTarget(validDepth, targetDepth)

dist    = abs(validDepth - targetDepth);
minDist = min(dist);
cand    = find(dist == minDist);

if numel(cand) == 1
    selRel = cand;
else
    [~, iMin] = min(validDepth(cand));   % tie-break: più superficiale
    selRel = cand(iMin);
end

depthSel = validDepth(selRel);

end