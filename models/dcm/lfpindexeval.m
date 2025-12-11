
function [total_elements,range_S,range_K] = lfpindexeval(data_trials)


num_elements_S  = length(data_trials(1).iLFP_S);
num_elements_K  = length(data_trials(1).iLFP_K);

total_elements  = num_elements_S + num_elements_K;
range_S         =1:num_elements_S;
range_K         =(num_elements_S + 1):total_elements;
return
if num_elements_S > 0
    range_S = sprintf('1:%d', num_elements_S);
else
    range_S = 'No Elements';
end

if num_elements_K > 0
    start_K = num_elements_S + 1;
    end_K = total_elements;
    range_K = sprintf('%d:%d', start_K, end_K);
else
    range_K = 'No Elements'; 
end