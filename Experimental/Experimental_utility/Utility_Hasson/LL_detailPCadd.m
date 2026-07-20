function [c_m, stats] = LL_detailPCadd(L_L_matrix, config_, tin1, tout1, tin2, tout2)

    [c_m, stats] = LL_detailPC(L_L_matrix, config_, tin1, tout1, tin2, tout2);

    stats = ensureStatFields(stats, 'S2K');
    stats = ensureStatFields(stats, 'K2S');
    stats = ensureStatFields(stats, 'DIFF');

    if isstruct(c_m)
        if isfield(c_m, 'S2K_vals')
            stats = addVectorStats(stats, c_m.S2K_vals, 'S2K');
        end
        if isfield(c_m, 'K2S_vals')
            stats = addVectorStats(stats, c_m.K2S_vals, 'K2S');
        end
        if isfield(c_m, 'S2K_vals') && isfield(c_m, 'K2S_vals')
            stats = addVectorStats(stats, c_m.K2S_vals - c_m.S2K_vals, 'DIFF');
        end
    end
end