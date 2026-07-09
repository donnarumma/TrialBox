% function Out_total = matrixArrange(In_session)

function Out_total = matrixArrange(In_session,par)

Out_total = struct();
nSession = numel(In_session);
for condIdx = 1:par.num_cond+1
    for dirIdx = 1:par.num_dir
        dirName = sprintf('dir%d', dirIdx);
        tempData = struct();
        for nSess = 1:nSession
            if isfield(In_session(nSess).cond(condIdx), dirName)
                tempData(nSess).(dirName) = In_session(nSess).cond(condIdx).(dirName);
            end
        end
        Out_total(condIdx).(dirName) = vertcat(tempData.(dirName));
    end
end
