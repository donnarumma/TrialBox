function out = mergeParams(def, user)
    out = def;
    fn = fieldnames(user);
    for i = 1:numel(fn)
        out.(fn{i}) = user.(fn{i});
    end
end
