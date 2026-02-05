function assert_valid_trial(trial, i)
% Validates data structure.
% Throws error if data don't satisfy the data contract.

    % --- required fields ---
    required_fields = {
        'trialTypeD'
        'trialTypeCon'
        'Spikes'
        'Manifold'
    };

    for f = 1:numel(required_fields)
        fname = required_fields{f};
        assert(isfield(trial, fname), ...
            'Trial %d: missing field "%s"', i, fname);
    end

    %  direction & condition -
    assert(~isempty(trial.trialTypeD), ...
        'Trial %d: trialTypeD is empty', i);

    assert(isscalar(trial.trialTypeD) && isnumeric(trial.trialTypeD), ...
        'Trial %d: trialTypeD must be numeric scalar', i);

    assert(~isempty(trial.trialTypeCon), ...
        'Trial %d: trialTypeCon is empty', i);

    assert(isscalar(trial.trialTypeCon) && isnumeric(trial.trialTypeCon), ...
        'Trial %d: trialTypeCon must be numeric scalar', i);

    % spike
    assert(~isempty(trial.Spikes), ...
        'Trial %d: Spikes is empty', i);

    assert(isnumeric(trial.Spikes), ...
        'Trial %d: Spikes must be numeric', i);

    %  manifold
    assert(~isempty(trial.Manifold), ...
        'Trial %d: Manifold is empty', i);

    assert(isnumeric(trial.Manifold), ...
        'Trial %d: Manifold must be numeric', i);

end
