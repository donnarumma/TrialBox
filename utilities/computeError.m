function   out=computeError(output, target)
% function out=computeError(output, target)
output = output(:); % consider all dimensions in the same variance
target = target(:);
Nsamples=size(output,1);

target_m       =mean(target);
target_rep     =repmat(target_m,Nsamples,1);
tm_diff        =abs(target-target_rep);
tar_var        =tm_diff.*tm_diff;
sum_tar_var    =sum(tar_var,1);

output_m       =mean(output);
output_rep     =repmat(output_m,Nsamples,1);
om_diff        =abs(output-output_rep);
out_var        =om_diff.*om_diff;
sum_out_var    =sum(out_var,1);

err_diff       =abs(target-output);
err_var        =err_diff .* err_diff;
sum_err_var    =sum(err_var,1);

STD_target     =sqrt(sum_tar_var/Nsamples);
STD_output     =sqrt(sum_out_var/Nsamples);

STD            =sqrt(sum_err_var/Nsamples);
SEM            =STD/sqrt(Nsamples);

ME             =sum_err_var/Nsamples; 
ZE             =sum_err_var/STD_target;
SQE            =(1/2)*sum_err_var;
MSQE           =(1/2)*sum_err_var/Nsamples;
RMSE           =sqrt(sum_err_var/Nsamples);
out.ME         =ME;
out.ZE         =ZE;
out.SQE        =SQE;
out.RMSE       =RMSE;
out.MSQE       =MSQE;
out.STD        =STD;
out.SEM        =SEM;
out.STD_target =STD_target;
out.STD_output =STD_output;
return;