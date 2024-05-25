function   h=plot_Montage(data_trials,par)
% function h=plot_Montage(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig    =figure;
    par.hfig=hfig;
else
    hfig    =par.hfig;
end 
InField         = par.InField;
col             = par.col;
hmShow          = par.hmShow;

% nChannels x nTimes x nTrials
Img_rec         = cat(3,data_trials.(InField)); 
numObservations = length(data_trials);
idx             = randperm(numObservations,hmShow);

montage (Img_rec(:,:,idx), 'BorderSize', [1,1], 'BackgroundColor', col);

if nargout>0
    h=hfig;
end
%% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end