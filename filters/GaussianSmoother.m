function   [data_trials,out] = GaussianSmoother(data_trials, par)
% function [data_trials,out] = GaussianSmoother(data_trials, par)
%
% create a field with Gaussian smoothing for each trial
%
% Gaussian kernel smoothing of data across time.
%
% INPUTS:
%
% yIn      - input data (yDim x T)
% kernSD   - standard deviation of Gaussian kernel, in msec
% stepSize - time between 2 consecutive datapoints in yIn, in msec
%
% OUTPUT:
%
% yOut     - smoothed version of yIn (yDim x T)
%
% OPTIONAL ARGUMENT:
%
% causal   - logical indicating whether temporal smoothing should
%            include only past data (true) or all data (false)
%
% @ from 2009 Byron Yu -- byronyu@stanford.edu
% 
% 2024 02 04: @CONAN - COgnitioN in ActioN
execinfo=par.exec;
if ~isempty(execinfo); t=tic; end

causal    = par.causal; 
kernSD    = par.kernSD;
stepSize  = par.stepSize;
nSD       = par.nSD;
InField   = par.InField;
OutField  = par.OutField;
Ntrials   = length(data_trials);

% Filter half length
% Go nSD standard deviations out
fltHL = ceil(nSD * kernSD / stepSize);

% Length of flt is 2*fltHL + 1
flt = normpdf(-fltHL*stepSize : stepSize : fltHL*stepSize, 0, kernSD);

if causal
    flt(1:fltHL) = 0;
end

for it=1:Ntrials 
    yIn       = data_trials(it).(InField);
    if (kernSD == 0) || (size(yIn, 2)==1)
        yOut = yIn;
    else
        [varDim, T] = size(yIn);
        yOut      = nan(varDim, T);
        
        % Normalize by sum of filter taps actually used
        nm = conv(flt, ones(1, T));
        
        for i = 1:varDim
            ys = conv(flt, yIn(i,:)) ./ nm;
            % Cut off edges so that result of convolution is same length 
            % as original data
            yOut(i,:) = ys(fltHL+1:end-fltHL);
        end
    end
    data_trials(it).(OutField)=yOut;
    data_trials(it).(['time' InField])=data_trials(it).(['time' OutField]);
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end
