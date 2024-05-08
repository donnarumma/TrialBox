function [o] = plotMeanAndVariancePlus(means,sigmas,params)
% function [] = plotMeanAndVariancePlus(means,sigmas,params)
% params.c1          = 'k';              % colour of the mean line
% params.c2          = [0.95,0.95,0];    % colour of the vertical lines
% params.c3          = 'r';              % colour of the sigmas lines
% params.hm          = 15;               % how many repetitions of vertical lines (increase if you see blank spaces)
% params.lw1         = 3;                % width of the vertical lines
% params.lw2         = 5;                % width of the mean line
% params.lw3         = 2;                % width of the sigmas line
% params.low         = -inf;             % inferior saturation
% params.high        = +inf;             % superior saturation
% params.fill        = 0;                % use fill (1 yes, 0 no)
try
    times = params.times;  % x values of the plot
catch
    times = 1:length(means);
end
try
    lw1 = params.lw1;
catch
    lw1 = 3;   % width of the vertical lines
end
try
    lw2 = params.lw2;
catch
    lw2         = 5;    % width of the mean line
end
try
    lw3 = params.lw3;
catch
    lw3 = 2;    % width of the sigmas line
end
try
    c1 = params.c1;
catch
    c1 ='k';              % colour of the mean line
end
try
    c2 = params.c2;
catch
    c2 = [0.95,0.95,0];    % colour of the vertical lines
end
try
    c3 = params.c3;
catch
    c3 = 'r';              % colour of the sigmas lines
end
try
    hm = params.hm;
catch
    hm          = 15;                   % how many repetitions (increase if you see blank spaces)
end
try
    low = params.low;
catch
    low = -inf;
end
try 
    high = params.high;
catch
    high = +inf;
end
try
    fill_option = params.fill;
catch
    fill_option = 0;
end
try
    hfig = params.hfig;
catch
    hfig =[];
end
if isempty(hfig)
    hfig=figure; hold on; box on; grid on;
    params.hfig=hfig;
end
T        =length(times)-1;
out.hbars=gobjects(T,hm);

lowlevel                 =means-sigmas;
lowlevel ( lowlevel<low) =low;
highlevel                =means+sigmas;
highlevel(highlevel>high)=high;

if fill_option
    out.bars = fill([times(:)' fliplr(times(:)')], [highlevel(:)' fliplr(lowlevel(:)')],c2);
    set(out.bars,'Edgecolor',c2);
else
    for i=1:T
        t   = linspace(    times(i),    times(i+1),hm);
        xd  = linspace( lowlevel(i), lowlevel(i+1),hm);
        xu  = linspace(highlevel(i),highlevel(i+1),hm);
        for k=1:hm
            out.hbars(i,k)=plot([t(k),t(k)],[xd(k),xu(k)],'color',c2,'LineWidth',lw1);
        end
    end
end
out.hmean   = plot(times,    means,'linewidth',lw2,'color',c1);
out.hinf    = plot(times, lowlevel,'linewidth',lw3,'color',c3);
out.hsup    = plot(times,highlevel,'linewidth',lw3,'color',c3);
out.hfig    = hfig;
if nargout > 0
    o=out;
end
return
        