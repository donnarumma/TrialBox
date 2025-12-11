function L = build_L_matrix(channel_source_map)
% build_L_matrix  Constructs a simple channel-to-source mapping matrix L
%
% FORMAT L = build_L_matrix(channel_source_map)
%
% channel_source_map : vector of length n_channels, with values in 1..n_sources
%
% OUTPUT:
% L : [n_channels x n_sources] lead field matrix
%     Each row maps one channel to one source

n_channels  = length(channel_source_map);
n_sources   = max(channel_source_map);

L = zeros(n_channels, n_sources);

for c = 1:n_channels
    s = channel_source_map(c);
    L(c, s) = 1;
end

end