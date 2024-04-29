function class_names = findclass(data,startClass)

trialnames          = {data.trialName};
[idx,idx_types]     = unique([data.trialType]);
idx_not             = setdiff(startClass,idx);
name = trialnames(idx_types);

class_names = cell(length(startClass),1);
for i=1:length(idx)
    class_names{idx(i)}    = name{i};
end
for j=1:length(idx_not)
    class_names{idx_not(j)} = 'NO';
end
