function [corr_avg] = f_corr(lag_s_X, lag_s_Y, varargin)
 
lag_s_X_3d=cell2mat(permute(lag_s_X,[3, 2, 1]))
lag_s_X_E= squeeze(mean(lag_s_X_3d, 2))';

lag_s_Y_3d=cell2mat(permute(lag_s_Y,[3, 2, 1]))
lag_s_Y_E= squeeze(mean(lag_s_Y_3d, 2))';

corr_avg=corr(lag_s_X_E',lag_s_Y_E')

end