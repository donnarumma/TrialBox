%function [PV_out,PV_session,lenTot,PV_outS,PV_outK,PV_name] = PVExtrat4eval(Move,lenM,par)

function [PV_out,lenTot,PV_outS,PV_outK,PV_name,PV_S,PV_K] = PVExtrat4eval(Move,lenM,par)


% parameters
num_cond    = par.num_cond;
num_dir     = par.num_dir;
idsess      = par.idsess;
fieldname   = par.fieldname;
time        = par.time;
windowSize  = par.windowSize;

MS              = struct();
MK              = struct();
TimePV    = struct();
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
        TimePV(cd).(strcat('dir',num2str(ndir)))  = timeSK;
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

        timeVel     = TimePV(cd).(strcat('dir', num2str(ndir)));
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
            MS_PV = sqrt(MS_vx.^2 + MS_vy.^2);


            ofilter         = 3; % Filter order Savitzky-Golay
            %% moving average filter
            MS_MA(k).Vel    = movmean(MS_PV, windowSize);

            if ~isnan(MS_PV)
                %% Savitzky-Golay filter
                MS_SG(k).Vel = sgolayfilt(MS_PV, ofilter, windowSize+1);
            else
                MS_SG(k).Vel = 0;
            end

            MK_dx            = diff(MK_data_x);
            MK_dy            = diff(MK_data_y);

            MK_vx            = MK_dx./dt;
            MK_vy            = MK_dy./dt;
            MK_PV      = sqrt(MK_vx.^2 + MK_vy.^2);

            %% moving average filter
            MK_MA(k).Vel     = movmean(MK_PV, windowSize);

            if ~isnan(MK_PV)
                %% Savitzky-Golay filter
                MK_SG(k).Vel = sgolayfilt(MK_PV, ofilter, windowSize+1);
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

MS_smoothSG_Max = struct();
MK_smoothSG_Max = struct();

for cd=1:num_cond
    for ndir = 1:num_dir
        time_app    = TimePV(cd).(strcat('dir', num2str(ndir)));
        Vs          = MS_smoothSG(cd).(strcat('dir',num2str(ndir)));
        Vk          = MK_smoothSG(cd).(strcat('dir',num2str(ndir)));
        Vs_max      = NaN(lenM(cd,ndir),1);
        Vk_max      = NaN(lenM(cd,ndir),1);
        for k=1:lenM(cd,ndir)
            TargOnset_ind   = find(time_app(k).alig>=0,1,'First');
            Vs_app          = Vs(k).Vel;
            Vk_app          = Vk(k).Vel;
            if Vs_app~=0
                Vs_max(k,1) = max(Vs_app(TargOnset_ind:end));
            else
                Vs_max(k,1) = 0;
            end
            if Vk_app~=0
                Vk_max(k,1) = max(Vk_app(TargOnset_ind:end));
            else
                Vk_max(k,1) = 0;
            end
        end
        MS_smoothSG_Max(cd).(strcat('dir',num2str(ndir))) = Vs_max;
        MK_smoothSG_Max(cd).(strcat('dir',num2str(ndir))) = Vk_max;
    end
end


PV_outS       = MS_smoothSG_Max;
PV_outK       = MK_smoothSG_Max;
PV_name       = M_name;
PV_S          = MS_smoothMA;
PV_K          = MK_smoothMA;

% Output
PV_out(1)     = MS_smoothSG_Max(1);
PV_out(2)     = MK_smoothSG_Max(2);
PV_out(3)     = MS_smoothSG_Max(3);
PV_out(4)     = MK_smoothSG_Max(3);

lenTot(idsess).len  = lenM;