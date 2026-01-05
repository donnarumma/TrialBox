function [A_sel,B_sel] = filter_data(A_struct, B_struct, config)
% describe

% INPUTS:
%   Trials       : Struct array (1 x N_trials).
%
% OUTPUT:
%   Filtered data for given subjects according to assigned condition and
%   directions

    dir_=config.dir;
    cond_= config.cond;
            % filter data according to direction and condition
    
    keep_A = false(1, numel(A_struct));
    keep_B = false(1, numel(B_struct));
    
    pippo_A = 0;
    pippo_B = 0;
    
    for i = 1:numel(A_struct)
        keep_A(i) = ...
          ismember(A_struct(i).trialTypeDir, dir_) && ...
          ismember(A_struct(i).trialTypeCond, cond_);
          pippo_A = pippo_A + 1;
        
    end
    
    for i = 1:numel(B_struct)
        keep_B(i) = ...
        ismember(B_struct(i).trialTypeDir, dir_) && ...
        ismember(B_struct(i).trialTypeCond, cond_);
        pippo_B = pippo_B + 1;
      
    end
    
    A_sel = A_struct(keep_A);
    B_sel = B_struct(keep_B);

end