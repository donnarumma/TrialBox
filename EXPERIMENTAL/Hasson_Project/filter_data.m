function [X_sel,Y_sel] = filter_data(X_struct, Y_struct, config)
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
    
    keep_X = false(1, numel(X_struct));
    keep_Y = false(1, numel(Y_struct));
    
    pippo_X = 0;
    pippo_Y = 0;
    
    for i = 1:numel(X_struct)
        keep_X(i) = ...
                ismember(X_struct(i).trialTypeDir, dir_) && ...
                ismember(X_struct(i).trialTypeCond, cond_);
          %ismember(X_struct(i).trialTypeDir, dir_) && ...
          %ismember(X_struct(i).trialTypeCond, cond_);
          pippo_X = pippo_X + 1;
        
    end
    
    for i = 1:numel(Y_struct)
        keep_Y(i) = ...
            any(ismember(Y_struct(i).trialTypeDir, dir_)) && ...
            any(ismember(Y_struct(i).trialTypeCond, cond_));
        pippo_Y = pippo_Y + 1;
      
    end
    
    X_sel = X_struct(keep_X);
    Y_sel = Y_struct(keep_Y);

end