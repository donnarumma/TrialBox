%function [Jxy_trajS,Jxy_trajK] = XYtrajectoryEval(MS,MK)


function [Jxy_trajS,Jxy_trajK] = XYtrajectoryEval(MS,MK)

num_cond = length(MS);
num_dir = length(fieldnames(MS));

traj_length_MS = struct();
traj_length_MK = struct();

for cd = 1:num_cond
    for ndir = 1:num_dir

        field = strcat('dir', num2str(ndir));

        % numero di trials per questa condizione/direzione
        num_trials = length(MS(cd).(field));

        for k = 1:num_trials

            % Estrazione dei dati per il trial k
            MS_data = MS(cd).(field)(k).XY;
            MK_data = MK(cd).(field)(k).XY;

            MS_data_x = MS_data(1,:);
            MS_data_y = MS_data(2,:);

            MK_data_x = MK_data(1,:);
            MK_data_y = MK_data(2,:);

            % Calcolo lunghezza traiettoria Monkey-S
            dx_MS = diff(MS_data_x);
            dy_MS = diff(MS_data_y);
            traj_length_MS(cd).(field)(k) = sum(sqrt(dx_MS.^2 + dy_MS.^2));

            % Calcolo lunghezza traiettoria Monkey-K
            dx_MK = diff(MK_data_x);
            dy_MK = diff(MK_data_y);
            traj_length_MK(cd).(field)(k) = sum(sqrt(dx_MK.^2 + dy_MK.^2));

        end
    end
end

Jxy_trajS = traj_length_MS;
Jxy_trajK = traj_length_MK;

end