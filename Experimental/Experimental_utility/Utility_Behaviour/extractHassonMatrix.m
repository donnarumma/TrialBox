% function extractHassonMatrix(H_matrix)

function [Cond1,Cond2,Cond3] = extractHassonMatrix(H_matrix,par)

numDirections = 8;
numConditions = 3;
Cond1 = NaN(1,numDirections);
Cond2 = NaN(1,numDirections);
Cond3 = NaN(1,numDirections);
for c = 1:numConditions
    for d = 1:numDirections
        % Trova riga con direzione d e condizione c
        idx = ([H_matrix.Direction] == d & [H_matrix.Condition] == c);

        m = H_matrix(idx).stats.(par.InField);
        
        if c == 1
            Cond1(1,d) = m;
        elseif c == 2
            Cond2(1,d) = m;
        else
            Cond3(1,d) = m;
        end
    end
end
