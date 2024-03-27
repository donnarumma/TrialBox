function   [h,out]=plot_trajectory3D(data_trials,par)
% function out=plot_trajectory3D(data_trials,par)
% keep : conditions to plot
% wd   : principal components to plot (must be three)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig    =figure;
    par.hfig=hfig;
else
    hfig    =par.hfig;
end

labels      = unique({data_trials.trialName});
all         = unique([data_trials.trialType]);

keep        = par.keep;
dimsToPlot  = par.wd;
explained   = par.explained;
cmaps       = par.cmaps;
center      = par.center;
istart      = par.istart;
iend        = par.iend;
axlab       = par.axlab;
InField     = par.InField;
nTrials     = length(data_trials);
    
lw          = 2.5;
ms          = 15;
        
hold on; box on; grid on;
set(gca,'fontsize',14)

pos = get(gcf, 'position');
set(hfig, 'position', [pos(1) pos(2) 1.3*pos(3) 1.3*pos(4)]);
mm  = 1;
try    
    v=zeros(nTrials,1);
    for ik=1:nTrials
      v(ik)=data_trials(ik).trialType;
    end
    mm=max(v);
end
  
for ik=mm
    out.h{ik}    =[];
    out.hp0{ik}  =[];
    out.hpc{ik}  =[];
    out.hpe{ik}  =[];
    out.PEnd{ik} =[];
end
for iTrial = 1:nTrials
    dat = data_trials(iTrial).(InField)(dimsToPlot,:);
    if isfield(data_trials(iTrial),'trialType')
        type=data_trials(iTrial).trialType;
        out.PEnd{type}=[out.PEnd{type};dat(1,iend), dat(2,iend), dat(3,iend)];
    end
end

for iTrial = 1:nTrials
    dat     = data_trials(iTrial).(InField)(dimsToPlot,:);
    type    = data_trials(iTrial).trialType;
    col     = cmaps(data_trials(iTrial).trialType,:);
    hp0=plot3(dat(1, istart), dat(2, istart), dat(3, istart),'o','markersize',ms,'MarkerEdgeColor','k','MarkerFaceColor',col);
    hpc=plot3(dat(1, center), dat(2, center), dat(3, center),'x','markersize',ms,'MarkerEdgeColor',col,'MarkerFaceColor',col,'linewidth',  5);
    hpe=plot3(dat(1,  iend ), dat(2,  iend ), dat(3,  iend ),'s','markersize',ms,'MarkerEdgeColor','k','MarkerFaceColor',col);
    hcur=plot3(dat(1,istart:iend), dat(2,istart:iend), dat(3,istart:iend), '.-', 'linewidth', lw, 'color', col);
    out.h  {type} =[out.h{type},  hcur];
    out.hp0{type} =[out.hp0{type},hp0 ];
    out.hpe{type} =[out.hpe{type},hpe ];
    out.hpc{type} =[out.hpc{type},hpc ];
    hold on;
end

axis square;

if isequal(InField, 'xorth')
    str1 = sprintf('$$\\tilde{\\mathbf %s}_{%d,t}$$', axlab, dimsToPlot(1));
    str2 = sprintf('$$\\tilde{\\mathbf %s}_{%d,t}$$', axlab, dimsToPlot(2));
    str3 = sprintf('$$\\tilde{\\mathbf %s}_{%d,t}$$', axlab, dimsToPlot(3));
else
    str1 = sprintf('$${\\mathbf %s}_{%d,:}$$', axlab, dimsToPlot(1));
    str2 = sprintf('$${\\mathbf %s}_{%d,:}$$', axlab, dimsToPlot(2));
    str3 = sprintf('$${\\mathbf %s}_{%d,:}$$', axlab, dimsToPlot(3));
end
xlabel(str1, 'interpreter', 'latex', 'fontsize', 24);
ylabel(str2, 'interpreter', 'latex', 'fontsize', 24);
zlabel(str3, 'interpreter', 'latex', 'fontsize', 24);

view(3); % view(70,16)
out.f=hfig;

ctd=all(~ismember(all,keep));
for ic=1:length(ctd)
    deleteCondition(out,ctd(ic));
end
labinds=[];
starts  =cell(size(labels));
ends    =cell(size(labels));
centers =cell(size(labels));
for iKeep=1:length(keep)
    ik           =keep(iKeep);
    labinds     =[labinds,   out.h{ik}];
    starts{ik}   =['start '  labels{ik}];
    ends{ik}     =['end '    labels{ik}];
    centers{ik}  =['center ' labels{ik}];
end
labinds=[labinds,[out.hp0{:}],[out.hpe{:}],[out.hpc{:}]];
labels =[labels(keep),starts,ends,centers];

legend(labinds,labels);
try
    ff=explained(dimsToPlot);
    d=sprintf('%3.0f',sum(ff));
    title(['Explained ' d '%']);
end
if nargout>0
    h=hfig;
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end