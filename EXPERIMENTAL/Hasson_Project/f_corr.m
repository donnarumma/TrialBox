function [corr_avg] = f_corr(lag_s_A, lag_s_B, varargin)
 
lag_s_A_3d=cell2mat(permute(lag_s_A,[3, 2, 1]))
lag_s_A_E= squeeze(mean(lag_s_A_3d, 2))';

lag_s_B_3d=cell2mat(permute(lag_s_B,[3, 2, 1]))
lag_s_B_E= squeeze(mean(lag_s_B_3d, 2))';

corr_avg=corr(lag_s_A_E',lag_s_B_E')



end