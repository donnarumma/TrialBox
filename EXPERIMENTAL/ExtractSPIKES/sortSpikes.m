% function data_sorted = sortSpikes(data)

function data_sorted = sortSpikes(data)

cond = [data.trialType2];   % 1,2,3
dir  = [data.trialType];    % 1..8

rep = zeros(size(dir));

for c = unique(cond)
    idxc = find(cond == c);
    rep(idxc) = ceil((1:numel(idxc))'/8);
end

[~, idx] = sortrows([cond(:) rep(:) dir(:)], [1 2 3]);

data_sorted = data(idx);
