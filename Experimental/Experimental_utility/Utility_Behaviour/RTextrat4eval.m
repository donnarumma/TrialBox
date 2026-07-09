%function [RT_out,RT_session,lenTot,RT_S,RT_K,RT_name] = RTextrat4eval(Move,lenM,par)

function [RT_out,lenTot,RT_S,RT_K,RT_name] = RTextrat4eval(Move,lenM,par)

% parameters
num_cond    = par.num_cond;
num_dir     = par.num_dir;
idsess      = par.idsess;
fieldname   = par.fieldname;

RT_S        = struct();
RT_K        = struct();
RT_NaN      = struct();
RT_name     = struct();
for cd=1:num_cond
    for ndir=1:num_dir
        Move_app                = Move(cd).(strcat('dir',num2str(ndir)));
        RT_Sdir                 = NaN(lenM(cd,ndir),1);
        RT_Kdir                 = NaN(lenM(cd,ndir),1);
        RT_trialID              = NaN(lenM(cd,ndir),1);
        for k = 1:lenM(cd,ndir)
            RT_Sdir(k,1)        = Move_app(k).(strcat(fieldname,'_S'));
            RT_Kdir(k,1)        = Move_app(k).(strcat(fieldname,'_K'));
            RT_trialID(k,1)     = Move_app(k).trialId;
        end
        RT_S(cd).(strcat('dir',num2str(ndir)))          = RT_Sdir;
        RT_K(cd).(strcat('dir',num2str(ndir)))          = RT_Kdir;
        RT_NaN(cd).(strcat('dir',num2str(ndir))).len    = length(RT_Sdir);
        RT_NaN(cd).(strcat('dir',num2str(ndir))).S      = [sum(isnan(RT_Sdir))];
        RT_NaN(cd).(strcat('dir',num2str(ndir))).K      = [sum(isnan(RT_Kdir))];
        RT_name(cd).(strcat('dir',num2str(ndir))).D     = ndir;
        RT_name(cd).(strcat('dir',num2str(ndir))).T     = RT_trialID;
    end
end

% Output
RT_out(1) = RT_S(1);
RT_out(2) = RT_K(2);
RT_out(3) = RT_S(3);
RT_out(4) = RT_K(3);

lenTot(idsess).len = lenM;