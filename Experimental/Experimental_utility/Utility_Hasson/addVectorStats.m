function stats = addVectorStats(stats, values, baseName)

    values = values(:);
    values = values(~isnan(values));

    stats.(['n_' baseName]) = numel(values);

    if isempty(values)
        stats.(['mean_' baseName]) = NaN;
        stats.(['std_'  baseName]) = NaN;
        stats.(['sem_'  baseName]) = NaN;
    else
        stats.(['mean_' baseName]) = mean(values, 'omitnan');
        stats.(['std_'  baseName]) = std(values, 0, 'omitnan');
        stats.(['sem_'  baseName]) = stats.(['std_' baseName]) / sqrt(numel(values));
    end
end