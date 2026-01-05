function [X_A,X_B, varargout] = ...
convert_to_lag_struct(A_sel,B_sel, config)

    % flag using neural data (hat-hat vel hat-obs)
    use_neural = config.use_neural
    % if use_neural
    %     disp('ciao')
    %    % assert(numel(varargin) == 2, ...
    %     %    'If neural data is used, pass exactly 2 extra inputs');
    % end

    if isempty(config.trial_length)
            trial_len = size(A_sel(1).Manifold, 2);
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
    n_trials = numel(A_sel);
    % filter data according to direction and condition
    for i=1:n_trials
        manif_A=A_sel(i).Manifold'
        manif_B=B_sel(i).Manifold'
        assert(size(manif_A,1) == trial_len && size(manif_B,1) == trial_len, ...
        'Manifold length mismatch with trial_length');
        if use_neural
           neural_A = A_sel(i).Spikes';
           neural_B = B_sel(i).Spikes';
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

            lag_manif_A{i} = cell(nB,1);
            lag_manif_B{i} = cell(nB,1);
        
            if use_neural
                lag_neural_A{i} = cell(nB,1);
                lag_neural_B{i} = cell(nB,1);
            end
  
        % add and concatenate
            for j = 1:nB
                s_b = blocks_idx(j, 1);
                e_b   = blocks_idx(j, 2);
                lag_manif_A{i}{j} = manif_A(s_b:e_b, :);
                lag_manif_B{i}{j} = manif_B(s_b:e_b, :);

                if use_neural
                    lag_neural_A{i}{j} = neural_A(s_b:e_b, :);
                    lag_neural_B{i}{j} = neural_B(s_b:e_b, :);
                end

            end
    end

    X_A=cell(nB,1)
    X_B=cell(nB,1)
    l_sb=sub_block_size
    stride_sb=sub_stride_abs

    for b = 1:nB
        X_A{b} = [];
        X_B{b} = [];
    
        for p = 1:n_trials
            emb_block_A = lag_manif_A{p}{b};   
            emb_block_B = lag_manif_B{p}{b};
            L_A = size(emb_block_A,1);   
            L_B = size(emb_block_B,1);
            for k_a=1:stride_sb:(L_A-l_sb+1)
                    idx_start_a=k_a;
                    idx_end_a=idx_start_a+l_sb-1;             
                    X_A{b} = [X_A{b}; mean(emb_block_A(idx_start_a:idx_end_a,:),1)];
            end
            for k_b=1:stride_sb:(L_B-l_sb+1)
                    idx_start_b=k_b;
                    idx_end_b=idx_start_b+l_sb-1;             
                    X_B{b} = [X_B{b}; mean(emb_block_B(idx_start_b:idx_end_b,:),1)];
             end
        end
    end
    if use_neural
        n_ch_A=size(neural_A,2)
        n_ch_B=size(neural_B,2)
        Y_A=cell(nB,n_ch_A)
        Y_B=cell(nB,n_ch_B)
        
        for b = 1:nB_ref
            for e =1:n_ch_A
                Y_A{b,e} = [];
                for p = 1:n_trials
                    current_block = lag_neural_A{p}{b}(:,e);   
                    L = length(current_block);
                    for  k=1:stride_sb:(L-l_sb+1)
                        idx_start=k;
                        idx_end=idx_start+l_sb-1;         
                        Y_A{b,e} = [Y_A{b,e}; mean(current_block(idx_start:idx_end))];
                    end
                end
            end
        end
    
        for b = 1:nB_ref
            for e =1:n_ch_B
                Y_B{b,e} = [];
                for p = 1:n_trials
                    current_block = lag_neural_B{p}{b}(:,e);   
                    L = length(current_block);
                    for  k=1:stride_sb:(L-l_sb+1)
                        idx_start=k;
                        idx_end=idx_start+l_sb-1;         
                        Y_B{b,e} = [Y_B{b,e}; mean(current_block(idx_start:idx_end))];
                    end
                end
            end
        end

    varargout{1} = Y_A;
    varargout{2} = Y_B;
end
end
