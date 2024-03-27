function   [h,out]=plot_trajectory2D(data_trials,par)
% function [h,out]=plot_trajectory2D(data_trials,par)
execinfo    = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig    = figure;
    par.hfig= hfig;
else
    hfig    = par.hfig;
end 
hold on; box on; grid on;
trialNames  = {data_trials.trialName};
trialTypes  = [data_trials.trialType];
[~,idx]     = unique(trialTypes);
labels      = trialNames(idx);
nClasses    = length(labels);
keep        = par.keep;
dimsToPlot  = par.wd;
explained   = par.explained;
cmaps       = par.cmaps;
center      = par.center;
istart      = par.istart;
iend        = par.iend;
axlab       = par.axlab;
InField     = par.InField;
    
lw          = 2.5;
ms          = 15;
        
nTrials     = length(data_trials);

set(gca,'fontsize',14)

pos         = get(gcf, 'position');
set(hfig, 'position', [pos(1) pos(2) 1.3*pos(3) 1.3*pos(4)]);
  
for iClass=1:nClasses
    out.h{iClass}    =[];
    out.hp0{iClass}  =[];
    out.hpc{iClass}  =[];
    out.hpe{iClass}  =[];
    out.PEnd{iClass} =[];
end

for iTrial = 1:nTrials
    dat     = data_trials(iTrial).(InField)(dimsToPlot,:);    
    type    = data_trials(iTrial).trialType;
    col     = cmaps(data_trials(iTrial).trialType,:);
    hp0     = plot(dat(1, istart), dat(2, istart), 'o','markersize',ms,'MarkerEdgeColor','k','MarkerFaceColor',col);
    hpc     = plot(dat(1, center), dat(2, center), 'x','markersize',ms,'MarkerEdgeColor',col,'MarkerFaceColor',col,'linewidth',  5);
    hpe     = plot(dat(1,  iend ), dat(2,  iend ), 's','markersize',ms,'MarkerEdgeColor','k','MarkerFaceColor',col);
    hcur    = plot(dat(1,istart:iend), dat(2,istart:iend), '.-', 'linewidth', lw, 'color', col);

    out.h  {type} =[out.h{type},  hcur];
    out.hp0{type} =[out.hp0{type},hp0 ];
    out.hpe{type} =[out.hpe{type},hpe ];
    out.hpc{type} =[out.hpc{type},hpc ];
end
axis square;

if isequal(InField, 'xorth')
    str1 = sprintf('$$\\tilde{\\mathbf %s}_{%d,t}$$', axlab, dimsToPlot(1));
    str2 = sprintf('$$\\tilde{\\mathbf %s}_{%d,t}$$', axlab, dimsToPlot(2));
else
    str1 = sprintf('$${\\mathbf %s}_{%d,:}$$', axlab, dimsToPlot(1));
    str2 = sprintf('$${\\mathbf %s}_{%d,:}$$', axlab, dimsToPlot(2));
end
xlabel(str1, 'interpreter', 'latex', 'fontsize', 24);
ylabel(str2, 'interpreter', 'latex', 'fontsize', 24);  
out.f=hfig;
%%
ctd=trialTypes(~ismember(trialTypes,keep));
for ic=1:length(ctd)
    deleteCondition(out,ctd(ic));
end
labinds =[]; hp0=[]; hpe=[]; hpc=[];
starts  =cell(size(labels));
ends    =cell(size(labels));
centers =cell(size(labels));
for iKeep=1:length(keep)
    ik           =keep(iKeep);
    labinds     =[labinds,   out.h{ik}(1)];
    hp0         =[hp0,     out.hp0{ik}(1)];
    hpe         =[hpe      out.hpe{ik}(1)];
    hpc         =[hpc      out.hpc{ik}(1)];
    starts{ik}  =['start '  labels{ik}];
    ends{ik}    =['end '    labels{ik}];
    centers{ik} =['center ' labels{ik}];
end
labinds=[labinds,hp0,hpe,hpc];
labels =[labels(keep),starts(keep),ends(keep),centers(keep)];
legend(labinds,labels,'Location','BestOutside');
try
    ff=explained(dimsToPlot);
    d=sprintf('%3.0f',sum(ff));
    title(['Explained ' d '%']);
end
if nargin>0
    h=hfig;
end
%% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end