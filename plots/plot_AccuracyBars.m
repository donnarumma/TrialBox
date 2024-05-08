function   h=plot_AccuracyBars(data_trials_prob_plot,par)
% function h=plot_AccuracyBars(data_trials_prob_plot,par)
%% plot ACCURACY with plot_EachDimBar
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig    =figure;
    par.hfig=hfig;
else
    hfig    =par.hfig;
end 
% meanData -> get class means
evaltime                            = par.evaltime;
SuccessField='success';
if isempty(evaltime)
    evaltime=data_trials_prob_plot(1).(['time' SuccessField]);
    evaltime=evaltime(end);
end
explained                           = par.explained;
channelSets                         = par.channelSets;
cmaps                               = par.cmaps;
nCols                               = par.nCols;
pms.meanData                        = meanDataParams;
pms.meanData.exec                   = [];
pms.meanData.trialTypes             = [data_trials_prob_plot.trialType];
pms.meanData.opt                    = [0,0,1]; % mean
pms.meanData.InField                = SuccessField;
pms.meanData.OutField               = 'accuracy';
class_acc_trials                    = meanData(data_trials_prob_plot,pms.meanData);
% meanData -> get all trial means
pms.meanData                        = meanDataParams;
pms.meanData.exec                   = [];
pms.meanData.opt                    = [0,0,1]; % mean
pms.meanData.trialTypes             = true(size([data_trials_prob_plot.trialType]));
pms.meanData.InField                = SuccessField;
pms.meanData.OutField               = 'accuracy';
acc_trials                          = meanData(data_trials_prob_plot,pms.meanData);   

% plot_EachDimBar
pms.plot_EachDimBar                 = plot_EachDimBarParams;
pms.plot_EachDimBar.exec            = [];
pms.plot_EachDimBar.hfig            = hfig;
pms.plot_EachDimBar.evaltime        = evaltime;
pms.plot_EachDimBar.novariance      = false;
pms.plot_EachDimBar.addbar          = false;
pms.plot_EachDimBar.cmaps           = cmaps;
pms.plot_EachDimBar.legplot         = 1;
pms.plot_EachDimBar.InField         = 'accuracy';
pms.plot_EachDimBar.novariance      = true;
pms.plot_EachDimBar.keep            = unique([class_acc_trials.trialType]);
pms.plot_EachDimBar.nCols           = nCols;
pms.plot_EachDimBar.ylabel          = '$acc$';
pms.plot_EachDimBar.YLIM            = [0-0.01,1+0.01];
pms.plot_EachDimBar.chanceline      = true;
% GRAPH SUBPLOT TITLES - one graph for each subset
nSets        = length(channelSets); %% number of graphs
str          = cell(nSets,1);
for iSet=1:nSets
    str{iSet}='';
end
for iSet=1:nSets
    channels     =channelSets{iSet};
    nChannels    =length(channels);
    % str{iSet}='';
    for ichannel=1:nChannels
        str{iSet} = sprintf('%s$${\\mathbf x}_{%d,t}$$',str{iSet},channels(ichannel));
        if ichannel<nChannels
             str{iSet} = sprintf('%s,',str{iSet});
        end 
    end
    if ~isempty(explained)
        str{iSet} = sprintf('%s (%2.1f)',str{iSet},sum(explained(channelSets{iSet})));
    end
end
pms.plot_EachDimBar.titles{1}   = 'ACCURACY';
pms.plot_EachDimBar.addbar      = acc_trials;%false;%data_trials_prob_sep(1);
hfig                            = plot_EachDimBar(class_acc_trials,pms.plot_EachDimBar);
sgtitle(['Classification at ' num2str(evaltime) ' s'])
if nargin > 0
    h=hfig;
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end