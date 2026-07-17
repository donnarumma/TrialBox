function stats = ensureStatFields(stats, baseName)

    if ~isfield(stats, ['std_' baseName])
        stats.(['std_' baseName]) = NaN;
    end
    if ~isfield(stats, ['sem_' baseName])
        stats.(['sem_' baseName]) = NaN;
    end
    if ~isfield(stats, ['n_' baseName])
        stats.(['n_' baseName]) = NaN;
    end
end