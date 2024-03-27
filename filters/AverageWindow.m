function   [data_trials,out] = AverageWindow(data_trials, par)
% function [data_trials,out] = AverageWindow(data_trials, par)
%
% seq = getSeq(dat, binWidth, ...)  
%
% Converts 0/1 spike trains into spike counts.
%
% INPUTS:
%
% data_Trials  structure whose nth entry (corresponding to the nth experimental
%               trial) has fields
%               (InField)-- matrix of the across all neurons.  Each row corresponds to a signal.
%                            Each column corresponds to 1 timestep (e.g. 1ms).
% binWidth    - bin width in msec
%
% OUTPUTS:
%
% data_trials  - data structure, whose nth entry (corresponding to
%               the nth experimental trial) has fields
%                 (OutField) with dimension (yDim x T) -- neural data
%
% OPTIONAL ARGUMENTS:
%
% useSqrt     - logical specifying whether or not to use square-root transform
%               on spike counts (default: true)
%
% @ from 2009 Byron Yu -- byronyu@stanford.edu
% 
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

useSqrt =par.useSqrt;
average =par.average;
binWidth=par.binWidth;
InField =par.InField;
OutField=par.OutField;


Ntrials =length(data_trials);
for itr = 1:Ntrials
    nChannels   = size(data_trials(itr).(InField), 1);
    if ~par.binWidth
        binWidth=size(data_trials(itr).(InField), 2);
    end
    nTimes      = floor(size(data_trials(itr).(InField), 2) / binWidth);
    
    inData       = data_trials(itr).(InField);
    data_trials(itr).(OutField) = nan(nChannels, nTimes);
    
    for iTime = 1:nTimes
        iStart = binWidth * (iTime-1) + 1;
        iEnd   = binWidth *  iTime;
        data_trials(itr).(OutField)(:,iTime) = sum(inData(:, iStart:iEnd), 2);
    end
    
    if useSqrt
        data_trials(itr).(OutField) = sqrt(data_trials(itr).(OutField));
    end
    if average
        data_trials(itr).(OutField) = data_trials(itr).(OutField)/binWidth;
    end
    time = data_trials(itr).(['time' InField]);
    deltat = (time(binWidth)-time(1))/2;
    data_trials(itr).(['time' OutField]) = linspace(time(1)+deltat,time(end)-deltat,nTimes);
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
end