function Data = build24BehaviourStruct(A, Data)
% build24BehaviorStruct
% Usage:
%   Data = build24BehaviorStruct(AE);
%   Data = build24BehaviorStruct(RT, Data);
%
% A must be a struct with fields:
%   A.S.Solo, A.S.Joint, A.K.Solo, A.K.Joint
%
% Output:
%   Data.S(1:24).<fieldname>
%   Data.K(1:24).<fieldname>
%
% Each call adds one field side-by-side inside the same 24 elements.

    if nargin < 1
        error('Usage: Data = build24BehaviorStruct(A) or Data = build24BehaviorStruct(A, Data)');
    end

    if ~isstruct(A) || ~isfield(A,'S') || ~isfield(A,'K')
        error('Input A must be a struct with fields A.S and A.K.');
    end
    if ~isfield(A.S,'Solo') || ~isfield(A.S,'Joint') || ~isfield(A.K,'Solo') || ~isfield(A.K,'Joint')
        error('A must contain A.S.Solo, A.S.Joint, A.K.Solo, and A.K.Joint.');
    end

    soloS  = A.S.Solo;
    jointS = A.S.Joint;
    soloK  = A.K.Solo;
    jointK = A.K.Joint;

    if numel(soloS) ~= 8 || numel(jointS) ~= 8 || numel(soloK) ~= 8 || numel(jointK) ~= 8
        error('Each input vector must have 8 elements.');
    end

    fn = inputname(1);
    if isempty(fn)
        fn = 'Value';
    end

    if nargin < 2 || isempty(Data)
        Data = struct();
        Data.S = repmat(struct(),1,24);
        Data.K = repmat(struct(),1,24);
    else
        if ~isfield(Data,'S') || ~isfield(Data,'K') || numel(Data.S) ~= 24 || numel(Data.K) ~= 24
            error('Existing Data must contain Data.S(1:24) and Data.K(1:24).');
        end
    end

    for i = 1:8
        Data.S(i).(fn) = soloS(i);
        Data.K(i).(fn) = NaN;
    end

    for i = 1:8
        Data.S(8+i).(fn) = NaN;
        Data.K(8+i).(fn) = soloK(i);
    end

    for i = 1:8
        Data.S(16+i).(fn) = jointS(i);
        Data.K(16+i).(fn) = jointK(i);
    end
end