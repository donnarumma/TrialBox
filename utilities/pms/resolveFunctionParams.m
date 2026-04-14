function par = resolveFunctionParams(functionName, parSource, varargin)
% resolveFunctionParams  Resolve default params, struct overrides, and name-value overrides.
%
% par = resolveFunctionParams(functionName)
% par = resolveFunctionParams(functionName, parSource)
% par = resolveFunctionParams(functionName, parSource, 'Name', value, ...)
%
% INPUT
%   functionName : char, string, or function handle. The corresponding
%                  params function is assumed to be [functionName 'Params'].
%   parSource    : optional struct with field overrides
%
% NOTES
%   Override precedence is:
%     defaults < parSource < name-value pairs

if nargin < 2 || isempty(parSource)
    parSource = struct();
elseif ~isstruct(parSource)
    varargin = [{parSource} varargin];
    parSource = struct();
end

paramsFunctionName = local_get_params_function_name(functionName);
paramsFunctionHandle = str2func(paramsFunctionName);
par = paramsFunctionHandle(parSource);

if isempty(varargin)
    return;
end

if mod(numel(varargin), 2) ~= 0
    error('%s expected name-value pairs after parSource.', mfilename);
end

for nameValueIndex = 1:2:numel(varargin)
    parameterName = varargin{nameValueIndex};
    parameterValue = varargin{nameValueIndex + 1};

    if isstring(parameterName)
        parameterName = char(parameterName);
    end

    if ~ischar(parameterName)
        error('%s expected parameter names to be char or string.', mfilename);
    end

    par.(parameterName) = parameterValue;
end
end


function paramsFunctionName = local_get_params_function_name(functionName)
if isa(functionName, 'function_handle')
    baseFunctionName = func2str(functionName);
elseif isstring(functionName)
    baseFunctionName = char(functionName);
elseif ischar(functionName)
    baseFunctionName = functionName;
else
    error('%s expected functionName to be char, string, or function handle.', mfilename);
end

if contains(baseFunctionName, '>')
    functionTokens = strsplit(baseFunctionName, '>');
    baseFunctionName = functionTokens{1};
end

[~, baseFunctionName, ~] = fileparts(baseFunctionName);

if endsWith(baseFunctionName, 'Params')
    paramsFunctionName = baseFunctionName;
else
    paramsFunctionName = [baseFunctionName 'Params'];
end
end
