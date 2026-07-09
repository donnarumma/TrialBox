% function [winner, p] = getBinom(wS, wK, total)

function [winner, p] = getBinom(wS, wK)
    total = wS + wK;
    if total == 0
        winner = "tie";
        p = NaN;
        return;
    end
    % Test binomiale a due code
    p = min(1, 2 * binocdf(min(wS, wK), total, 0.5));
    if wS > wK
        winner = "S";
    elseif wK > wS
        winner = "K";
    else
        winner = "tie";
    end
end
