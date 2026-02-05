function [X_,Y_, varargout] = ...
convert_to_lag_struct(X_sel,Y_sel, config)

    % flag using neural data (hat-hat vel hat-obs)
    use_neural = config.use_neural
    % if use_neural
    %     disp('ciao')
    %    % assert(numel(varargin) == 2, ...
    %     %    'If neural data is used, pass exactly 2 extra inputs');
    % end

    if isempty(config.trial_length)
            trial_len = size(X_sel(1).Manifold, 2);
    else
        trial_len=config.trial_length;
    end

     
    % Block parameters
    block_size=config.block_size
    block_stride_perc=config.block_stride
    % consideriamo di gestire casi patologici tipo valore non intero
    stride_abs   =floor(block_stride_perc * block_size);
    assert(stride_abs >= 1, 'stride_abs became 0: increase block_stride or block_size.');
    assert(block_size < trial_len, 'block_size must be smaller than trial_length');
    assert(block_stride_perc > 0 && block_stride_perc <= 1, 'block_stride must be in (0,1]');


  
    
    % "trial expolding" parameters
    sub_block_size=config.sub_block_size
    sub_block_stride_perc=config.sub_block_stride
    assert(sub_block_size <= block_size, 'sub_block_size must be <= block_size');
    sub_stride_abs   = floor(sub_block_stride_perc * sub_block_size);
    assert(sub_stride_abs >= 1,'sub_block_stride too small: stride became 0.');
    
    n_blocks = floor((trial_len - block_size) / stride_abs) + 1;
    
    % Pre-allocate
    % lag_manif_A = cell(n_trials,1);
    % lag_manif_B = cell( n_trials,1);
    % if use_neural
    %     lag_neural_A = cell(n_trials, 1);
    %     lag_neural_B = cell(n_trials, 1);
    % end
    n_trials = numel(X_sel);
    % filter data according to direction and condition
    for i=1:n_trials
        manif_X=X_sel(i).Manifold'
        manif_Y=Y_sel(i).Manifold'
        assert(size(manif_X,1) == trial_len && size(manif_Y,1) == trial_len, ...
        'Manifold length mismatch with trial_length');
        if use_neural
           neural_X = X_sel(i).Spikes';
           neural_Y = Y_sel(i).Spikes';
        end
        start_idx   = 1;
            %end_idx=trial_len_res;
        blocks_idx = [];
        while true
            end_idx = start_idx + block_size - 1;
            if end_idx > trial_len
               break
            end
                blocks_idx = [blocks_idx; start_idx end_idx];
                start_idx  = start_idx + stride_abs;
            end
     
            nB = size(blocks_idx, 1);
            if i == 1
                nB_ref = nB;
            else
                assert(nB == nB_ref,'Inconsistent number of blocks across trials');
            end

            lag_manif_X{i} = cell(nB,1);
            lag_manif_Y{i} = cell(nB,1);
        
            if use_neural
                lag_neural_X{i} = cell(nB,1);
                lag_neural_Y{i} = cell(nB,1);
            end
  
        % add and concatenate
            for j = 1:nB
                s_b = blocks_idx(j, 1);
                e_b   = blocks_idx(j, 2);
                lag_manif_X{i}{j} = manif_X(s_b:e_b, :);
                lag_manif_Y{i}{j} = manif_Y(s_b:e_b, :);

                if use_neural
                    lag_neural_X{i}{j} = neural_X(s_b:e_b, :);
                    lag_neural_Y{i}{j} = neural_Y(s_b:e_b, :);
                end

            end
    end

    X_=cell(nB,1)
    Y_=cell(nB,1)
    l_sb=sub_block_size
    stride_sb=sub_stride_abs

    for b = 1:nB
        X_{b} = [];
        Y_{b} = [];
    
        for p = 1:n_trials
            emb_block_X = lag_manif_X{p}{b};   
            emb_block_Y = lag_manif_Y{p}{b};
            L_X = size(emb_block_X,1);   
            L_Y = size(emb_block_Y,1);
            for k_a=1:stride_sb:(L_X-l_sb+1)
                    idx_start_a=k_a;
                    idx_end_a=idx_start_a+l_sb-1;             
                    X_{b} = [X_{b}; mean(emb_block_X(idx_start_a:idx_end_a,:),1)];
            end
            for k_b=1:stride_sb:(L_Y-l_sb+1)
                    idx_start_b=k_b;
                    idx_end_b=idx_start_b+l_sb-1;             
                    Y_{b} = [Y_{b}; mean(emb_block_Y(idx_start_b:idx_end_b,:),1)];
             end
        end
    end
    if use_neural
        n_ch_X=size(neural_X,2)
        n_ch_Y=size(neural_Y,2)
        Y_x=cell(nB,n_ch_X)
        Y_y=cell(nB,n_ch_Y)
        
        for b = 1:nB_ref
            for e =1:n_ch_X
                Y_x{b,e} = [];
                for p = 1:n_trials
                    current_block = lag_neural_X{p}{b}(:,e);   
                    L = length(current_block);
                    for  k=1:stride_sb:(L-l_sb+1)
                        idx_start=k;
                        idx_end=idx_start+l_sb-1;         
                        Y_x{b,e} = [Y_x{b,e}; mean(current_block(idx_start:idx_end))];
                    end
                end
            end
        end
    
        for b = 1:nB_ref
            for e =1:n_ch_Y
                Y_y{b,e} = [];
                for p = 1:n_trials
                    current_block = lag_neural_Y{p}{b}(:,e);   
                    L = length(current_block);
                    for  k=1:stride_sb:(L-l_sb+1)
                        idx_start=k;
                        idx_end=idx_start+l_sb-1;         
                        Y_y{b,e} = [Y_y{b,e}; mean(current_block(idx_start:idx_end))];
                    end
                end
            end
        end

    varargout{1} = Y_x;
    varargout{2} = Y_y;
end
end
