%function [JxyS_out,JxyK_out,Jxy_name,Jxy_time,RT_Jxy] = TrajXYextrat4eval(Move,par)

function [JxyS_out,JxyK_out,Jxy_name,Jxy_time,RT_Jxy] = TrajXYextrat4eval(Move,par)

% parameters
num_cond    = par.num_cond;
num_dir     = par.num_dir;
fieldname   = par.fieldname;
time        = par.time;

MS = struct();
MK = struct();
RT = struct();
TimeVelocity = struct();
M_name = struct();
% K = struct();
% K_trialID = struct();

for cd = 1:num_cond
    for ndir = 1:num_dir

        M_app = Move(cd).(strcat('dir',num2str(ndir)));

        % numero di trials disponibili per questa condizione/direzione
        num_trials = length(M_app);

        for k = 1:num_trials

            time_app = M_app(k).(time);

            MS_appX = double(M_app(k).(fieldname)(1,:));
            MS_appY = double(M_app(k).(fieldname)(2,:));

            MK_appX = double(M_app(k).(fieldname)(5,:));
            MK_appY = double(M_app(k).(fieldname)(6,:));

            % Salvataggio dati soggetto S
            MS(cd).(strcat('dir',num2str(ndir)))(k).XY = [MS_appX; MS_appY];

            % Salvataggio dati soggetto K
            MK(cd).(strcat('dir',num2str(ndir)))(k).XY = [MK_appX; MK_appY];

            % Tempo
            TimeVelocity(cd).(strcat('dir',num2str(ndir)))(k).time = time_app;

            % Informazioni sul trial
            M_name(cd).(strcat('dir',num2str(ndir)))(k).D = ndir;
            M_name(cd).(strcat('dir',num2str(ndir)))(k).T = M_app(k).trialId;

            % RT
            RT(cd).(strcat('dir',num2str(ndir)))(k).S = M_app(k).RT_S;
            RT(cd).(strcat('dir',num2str(ndir)))(k).K = M_app(k).RT_K;

        end
    end
end

JxyS_out    = MS;
JxyK_out    = MK;
Jxy_name    = M_name;
Jxy_time    = TimeVelocity;
RT_Jxy      = RT;

end

