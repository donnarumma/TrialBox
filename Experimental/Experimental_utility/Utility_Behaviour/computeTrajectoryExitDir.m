
function Result = computeTrajectoryExitDir(Jxy_session,par)

r_dim = par.r_dim;
n_sessions = length(Jxy_session);
monkey = {'S', 'K'};
Result = struct();

for s = 1:n_sessions
    Result(s).session = Jxy_session(s).session;
    for sj = 1:length(monkey)
        monkey_name = monkey{sj};
        for cond = 1:3
            for n_dir = 1:8
                dirName = sprintf('dir%d', n_dir);

                try
                    currtrials  = Jxy_session(s).(monkey_name)(cond).(dirName);
                    currTimes   = Jxy_session(s).time(cond).(dirName);
                    currTrId    = Jxy_session(s).name(cond).(dirName);
                    currRT      = Jxy_session(s).RT_Jxy(cond).(dirName);
                catch
                    continue;
                end

                n_trials = length(currtrials);
                trial_stats = struct();

                for t = 1:n_trials
                    data = currtrials(t).XY;
                    t_vec = currTimes(t).time;

                    idx_post_zero = find(t_vec >= 0);

                    info = struct();
                    info.target_dir = n_dir;

                    if isempty(idx_post_zero)
                        info.status = 'no_data_post_onset';
                        info.real_dir = NaN;
                        info.exit_time = NaN;
                        info.is_correct = false;
                    else
                        x = data(1, idx_post_zero);
                        y = data(2, idx_post_zero);
                        t_filtered = t_vec(idx_post_zero);

                        distances = sqrt(x.^2 + y.^2);
                        idx_exit = find(distances > r_dim, 1, 'first');

                        if isempty(idx_exit)
                            info.real_dir = NaN;
                            info.exit_time = NaN;
                            info.is_correct = false;
                            info.status = 'not_exit';
                        else
                            x_out = x(idx_exit);
                            y_out = y(idx_exit);
                            info.exit_time = t_filtered(idx_exit);

                            angle_rad = atan2(y_out, x_out);
                            angle_relative = pi/2 - angle_rad; % Ruota e inverte in senso orario

                            angle_real_deg = rad2deg(angle_relative);
                            angle_target_deg = (n_dir - 1) * 45;

                            info.angular_error = mod((angle_real_deg - angle_target_deg) + 180, 360) - 180;
                            % sector1 = mod(round(-(angle_rad - pi/2) / (pi/4)), 8);

                            sector = mod(floor((angle_relative + pi/8) / (pi/4)), 8);

                            if sector == 0, sector = 8; end

                            final_sector = sector + 1;
                            if final_sector > 8, final_sector = 1; end

                            info.real_dir = final_sector;
                            info.is_correct = (info.real_dir == n_dir);
                            info.status = 'exit';
                        end
                        info.TrialID = currTrId(t).T;
                        info.CMT      = info.exit_time - currRT(t).(monkey_name); % Central Movement Time
                    end
                    trial_stats(t).info = info;
                end
                Result(s).(monkey_name)(cond).(dirName) = trial_stats;
            end
        end
    end
end