function   h=plot_ellipsoid2D(seq,selected,indsT,cmaps,par)
% function h=plot_ellipsoid2D(seq,selected,indsT,cmaps,par)
% 2D basis
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig    =figure;
    par.hfig=hfig;
else
    hfig    =par.hfig;
end 
xspec = par.InField;

try
    xtkl    = seq.(['time' (xspec)]);
catch
    xtkl    = getXTKL(par);
end
pc      =[1,0; ...
          0,1];
hold on; box on; grid on; %axis('equal');
axlab='x';
str1 = sprintf('$$\\tilde{\\mathbf %s}_{%d,t}$$', axlab, 1);
str2 = sprintf('$$\\tilde{\\mathbf %s}_{%d,t}$$', axlab, 2);
Ncond           =length(selected);
try
    explained=par.explained;
    str1 = sprintf('%s (%2.1f)',str1,explained(1));
    str2 = sprintf('%s (%2.1f)',str2,explained(2));
catch
end
try
    dimsToPlot = par.wd;
catch
    dimsToPlot = 1:par.xDim;
end

Tot     =length(indsT);
he=gobjects(Tot,Ncond);
h =gobjects(Tot,Ncond);

%% plot ellipse
all_labels=cell(Tot,Ncond);
for it=1:length(indsT)
    indT=indsT(it);
    for isel=1:length(selected)
        sel                 = selected(isel);
        [X0,~,labnames]     = getPointsPCA(seq,sel,indT,dimsToPlot,xspec);
        all_labels(it,isel) = unique(labnames);
        [mm,ss,~]=CI_compute(X0,par);
  
        xy=mean(X0); timefield   = linspace(0,2*pi,30); data = pc*[ss(1)*cos(timefield);ss(2)*sin(timefield)];
        try
            he(it,isel)=plot3(xtkl(indT)*ones(size(data(1,:))),xy(1)+data(1,:),xy(2)+data(2,:),'color',cmaps(sel,:),'linewidth',3);
        catch
        end
        h(it,isel)=plot3(xtkl(indT)*ones(size(mm(:,1))),mm(1),mm(2),'o','color',cmaps(sel,:),'linewidth',3);
    end
end
xlabel('t (ms)','fontsize', 24);
ylabel(str1, 'interpreter', 'latex', 'fontsize', 24);
zlabel(str2, 'interpreter', 'latex', 'fontsize', 24);

p=[0,0,1000,600];
fs=16;
set(gca,'fontsize',fs);
set(hfig, 'PaperPositionMode','auto');
set(hfig,'Position',p);

legend(h(1,:),all_labels(1,:),'Location','Best');
view(-10,40)
if nargout>0
    h=hfig;
end
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end