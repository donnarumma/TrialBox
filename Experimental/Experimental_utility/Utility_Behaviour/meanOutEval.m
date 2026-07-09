% function [Out_Mean,Out_Std,Out_SE,max_y] = meanOutEval(In_session)

function [Out_Mean,Out_Std,Out_SE,max_y] = meanOutEval(In_session,par)

In_Mean = struct();
In_Std = struct();
In_SE = struct();
for icond = 1:par.num_cond+1
    for ndir = 1:par.num_dir
        In_calc = In_session(icond).(strcat('dir',num2str(ndir)));
        In_Mean(icond).(strcat('dir',num2str(ndir))) = mean(In_calc);
        In_Std(icond).(strcat('dir',num2str(ndir))) = std(In_calc);
        In_SE(icond).(strcat('dir',num2str(ndir))) = std(In_calc)/length(In_calc);
    end
end
Out_Mean    = (cell2mat(struct2cell(In_Mean')))';
Out_Std     = (cell2mat(struct2cell(In_Std')))';
Out_SE      = (cell2mat(struct2cell(In_SE')))';
max_y       = max(max(Out_Mean+2*Out_Std));
