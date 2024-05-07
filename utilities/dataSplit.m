function   [data, out] = dataSplit(data,params)
% function [data, out] = dataSplit(data,params)

execinfo=params.exec;
if ~isempty(execinfo); t=tic; end
TrainPercentage     = params.TrainPercentage;
TrainPercentageV    = params.TrainPercentageV;

labels              = [data.trialType]';

% Split_train_test
NTrials             = length(data);
TrSize              = floor((NTrials*TrainPercentage)/(2*100))*2;
TsSize              = NTrials-TrSize;
part_train_test     = cvpartition(labels,'HoldOut',TsSize,'Stratify',true);

% log_train           = training(data_partition,1);
log_train           = training(part_train_test,1);
log_test            = test(part_train_test, 1);
labels_train_val    = labels(log_train);

log_valid           = false(size(log_train));
part_train_valid    = []; 
if TrainPercentageV>0
    % Split_train_validation
    NValTrials          = sum(log_train);
    TrVSize             = floor((NValTrials*TrainPercentageV)/(2*100))*2;
    ValSize             = NValTrials-TrVSize;

    part_train_valid    = cvpartition(labels_train_val,'HoldOut',ValSize,'Stratify',true);
    
    logon_idx_valid     = test(part_train_valid, 1);
    idx_train           = find(log_train);
    for idt=1:length(idx_train)
        if logon_idx_valid(idt)
            id           = idx_train(idt);
            log_valid(id)=true;
        end
    end
end

for it=1:NTrials
    data(it).train  = log_train(it);
    data(it).test   = log_test(it);
    data(it).valid  = log_valid(it);
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end
out.part_train_test     = part_train_test;
out.part_train_valid    = part_train_valid;