%function [PVT_out,Time_session,lenTot,PVT_outS,PVT_outK,Time_name] = PVTExtrat4eval(Move,lenM,par)

function [PVT_out,lenTot,PVT_outS,PVT_outK,Time_name] = PVTExtrat4eval(Move,lenM,par)


% parameters
num_cond        = par.num_cond;
num_dir         = par.num_dir;
idsess          = par.idsess;
fieldname       = par.fieldname;
time            = par.time;
windowSize      = par.windowSize;

MS              = struct();
MK              = struct();
PVT             = struct();
M_name          = struct();
for cd=1:num_cond
    for ndir=1:num_dir
        M_app   = Move(cd).(strcat('dir',num2str(ndir)));
        MS_app  = struct();
        MK_app  = struct();
        timeSK  = struct();
        trialID = NaN(lenM(cd,ndir));
        for k = 1:lenM(cd,ndir)
            time_app        = M_app(k).(time);
            MS_appX         = double(M_app(k).(fieldname)(1,:));
            MS_appY         = double(M_app(k).(fieldname)(2,:));
            MK_appX         = double(M_app(k).(fieldname)(5,:));
            MK_appY         = double(M_app(k).(fieldname)(6,:));
            MS_app(k).XY    = [MS_appX;MS_appY];
            MK_app(k).XY    = [MK_appX;MK_appY];
            timeSK(k).alig  = time_app;
            trialID(k,1)    = M_app(k).trialId;
        end
        MS(cd).(strcat('dir',num2str(ndir)))            = MS_app;
        MK(cd).(strcat('dir',num2str(ndir)))            = MK_app;
        PVT(cd).(strcat('dir',num2str(ndir)))  = timeSK;
        M_name(cd).(strcat('dir',num2str(ndir))).D      = ndir;
        M_name(cd).(strcat('dir',num2str(ndir))).T      = trialID;
        % K(cd).(strcat('dir',num2str(ndir)))           = k;
        % K_trialID(cd).(strcat('dir',num2str(ndir)))   = trialK;
    end
end

MS_smoothMA = struct();
MS_smoothSG = struct();
MK_smoothMA = struct();
MK_smoothSG = struct();
for cd=1:num_cond
    for ndir = 1:num_dir
        MS_data     = MS(cd).(strcat('dir', num2str(ndir)));

        MK_data     = MK(cd).(strcat('dir', num2str(ndir)));

        timeVel     = PVT(cd).(strcat('dir', num2str(ndir)));
        MS_MA       = struct();
        MS_SG       = struct();
        MK_MA       = struct();
        MK_SG       = struct();
        for k=1:lenM(cd,ndir)
            MS_data_x   = MS_data(k).XY(1,:);
            MS_data_y   = MS_data(k).XY(2,:);
            MK_data_x   = MK_data(k).XY(1,:);
            MK_data_y   = MK_data(k).XY(2,:);
    
            dt          = diff(timeVel(k).alig);

            MS_dx       = diff(MS_data_x);
            MS_dy       = diff(MS_data_y);

            MS_vx       = MS_dx./dt;
            MS_vy       = MS_dy./dt;
            MS_velocity = sqrt(MS_vx.^2 + MS_vy.^2);


            ofilter         = 3; % Filter order Savitzky-Golay
            %% moving average filter
            MS_MA(k).Vel    = movmean(MS_velocity, windowSize);

            if ~isnan(MS_velocity)
                %% Savitzky-Golay filter
                MS_SG(k).Vel = sgolayfilt(MS_velocity, ofilter, windowSize+1);
            else
                MS_SG(k).Vel = 0;
            end

            MK_dx            = diff(MK_data_x);
            MK_dy            = diff(MK_data_y);

            MK_vx            = MK_dx./dt;
            MK_vy            = MK_dy./dt;
            MK_velocity      = sqrt(MK_vx.^2 + MK_vy.^2);

            %% moving average filter
            MK_MA(k).Vel     = movmean(MK_velocity, windowSize);

            if ~isnan(MK_velocity)
                %% Savitzky-Golay filter
                MK_SG(k).Vel = sgolayfilt(MK_velocity, ofilter, windowSize+1);
            else
                MK_SG(k).Vel = 0;
            end
            MS_smoothMA(cd).(strcat('dir',num2str(ndir))) = MS_MA;
            MS_smoothSG(cd).(strcat('dir',num2str(ndir))) = MS_SG;
            MK_smoothMA(cd).(strcat('dir',num2str(ndir))) = MK_MA;
            MK_smoothSG(cd).(strcat('dir',num2str(ndir))) = MK_SG;
        end
    end
end

TS_smoothSG_Max = struct();
TK_smoothSG_Max = struct();

for cd=1:num_cond
    for ndir = 1:num_dir
        time_app    = PVT(cd).(strcat('dir', num2str(ndir)));
        Vs          = MS_smoothSG(cd).(strcat('dir',num2str(ndir)));
        Vk          = MK_smoothSG(cd).(strcat('dir',num2str(ndir)));
        Vs_max      = NaN(lenM(cd,ndir),1);
        Vk_max      = NaN(lenM(cd,ndir),1);
        Ts_max      = NaN(lenM(cd,ndir),1);
        Tk_max      = NaN(lenM(cd,ndir),1);
        for k=1:lenM(cd,ndir)
            time_appK = time_app(k).alig;
            TargOnset_ind   = find(time_app(k).alig>=0,1,'First');
            time_TgOn = time_appK(TargOnset_ind:end);
            Vs_app          = Vs(k).Vel;
            Vk_app          = Vk(k).Vel;
            if Vs_app~=0
                [Vs_max(k,1),Ts_maxind] = max(Vs_app(TargOnset_ind:end));
                Ts_max(k,1)             = time_TgOn(Ts_maxind);
            else
                Vs_max(k,1) = 0;
                Ts_maxind   = 0;
                Ts_max(k,1) = 0;
            end
            if Vk_app~=0
                [Vk_max(k,1),Tk_maxind]     = max(Vk_app(TargOnset_ind:end));
                Tk_max(k,1)                 = time_TgOn(Tk_maxind);
            else
                Vk_max(k,1) = 0;
                Tk_maxind   = 0;
                Tk_max(k,1) = 0;
            end
        end
        TS_smoothSG_Max(cd).(strcat('dir',num2str(ndir))) = Ts_max;
        TK_smoothSG_Max(cd).(strcat('dir',num2str(ndir))) = Tk_max;
    end
end

PVT_outS               = TS_smoothSG_Max;
PVT_outK               = TK_smoothSG_Max;
Time_name               = M_name;
% Output
PVT_out(1)             = TS_smoothSG_Max(1);
PVT_out(2)             = TK_smoothSG_Max(2);
PVT_out(3)             = TS_smoothSG_Max(3);
PVT_out(4)             = TK_smoothSG_Max(3);

lenTot(idsess).len      = lenM;