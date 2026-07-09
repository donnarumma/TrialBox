%function [ST_out,lenTot,ST_S,ST_K,ST_name] = STextrat4eval(Move,lenM,par)

function [ST_out,lenTot,ST_S,ST_K,ST_name] = STextrat4eval(Move,lenM,par)

% parameters
num_cond = par.num_cond;
num_dir = par.num_dir;
idsess = par.idsess;
fieldname = par.fieldname;

ST_S = struct();
ST_K = struct();
ST_NaN = struct();
ST_name = struct();
for cd=1:num_cond
    for ndir=1:num_dir
        Move_app = Move(cd).(strcat('dir',num2str(ndir)));
        ST_Sdir = NaN(lenM(cd,ndir),1);
        ST_Kdir = NaN(lenM(cd,ndir),1);
        ST_trialID = NaN(lenM(cd,ndir),1);
        for k = 1:lenM(cd,ndir)
            ST_Sdir(k,1) = Move_app(k).(strcat(fieldname));
            ST_Kdir(k,1) = Move_app(k).(strcat(fieldname));
            ST_trialID(k,1)   = Move_app(k).trialId;
        end
        ST_S(cd).(strcat('dir',num2str(ndir))) = ST_Sdir;
        ST_K(cd).(strcat('dir',num2str(ndir))) = ST_Kdir;
        ST_NaN(cd).(strcat('dir',num2str(ndir))).len = length(ST_Sdir);
        ST_NaN(cd).(strcat('dir',num2str(ndir))).S = [sum(isnan(ST_Sdir))];
        ST_NaN(cd).(strcat('dir',num2str(ndir))).K = [sum(isnan(ST_Kdir))];
        ST_name(cd).(strcat('dir',num2str(ndir))).D = ndir;
        ST_name(cd).(strcat('dir',num2str(ndir))).T = ST_trialID;
    end
end

% Output
ST_out(1) = ST_S(1);
ST_out(2) = ST_K(2);
ST_out(3) = ST_S(3);
ST_out(4) = ST_K(3);

lenTot(idsess).len = lenM;