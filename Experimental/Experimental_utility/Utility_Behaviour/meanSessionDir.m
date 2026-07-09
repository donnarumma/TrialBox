% function [data_outMEAN,data_outSTD] = meanSessionDir(data_session,par)

function [data_outMEAN,data_outSTD] = meanSessionDir(data_session,par)

data_outMEAN = struct();
data_outSTD = struct();
session_name = par.SessName;
for n_sess = 1:length(data_session)
    data_outMEAN(n_sess).name = session_name{n_sess};
    data_outSTD(n_sess).name = session_name{n_sess};
    data_appSess = data_session(n_sess);
    for n_cond = 1:par.num_cond+1
        for dirIdx = 1:par.num_dir
            if par.millisec
                dirName = sprintf('dir%d', dirIdx);
                data_outMEAN(n_sess).cond(n_cond).(dirName) = 1000*mean(data_appSess.cond(n_cond).(dirName));
                data_outSTD(n_sess).cond(n_cond).(dirName) = 1000*(std(data_appSess.cond(n_cond).(dirName))/length(data_appSess.cond(n_cond).(dirName)));
            else 
                dirName = sprintf('dir%d', dirIdx);
                data_outMEAN(n_sess).cond(n_cond).(dirName) = mean(data_appSess.cond(n_cond).(dirName));
                data_outSTD(n_sess).cond(n_cond).(dirName) = std(data_appSess.cond(n_cond).(dirName))/length(data_appSess.cond(n_cond).(dirName));
            end
        end
    end
end
