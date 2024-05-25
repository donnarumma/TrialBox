function   h = plot_Atoms(data_trials,par)
% function h = plot_Atoms(data_trials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end
if isempty(par.hfig) || ~isvalid(par.hfig)
    hfig    =figure;
    par.hfig=hfig;
else
    hfig    =par.hfig;
end 
Aname               = par.Aname;
% hold on; box on; grid on;
InField             = par.InField;
nRows               = par.nRows;
xfld                = par.xfld;
times               = data_trials(1).([xfld, InField]);
[nChannels, nTimes] = size(data_trials(1).(InField));
nAtoms              = length(data_trials);
ylines              = 1:nChannels:nChannels*(nAtoms+1);
ylinesup            = ylines(2:end);
for iAtom=1:nAtoms
    labnames{iAtom} = sprintf('%s%g',Aname,iAtom);
end
cmaps               = linspecer(nAtoms);

% nChannels x nTiems x nAtoms
ID = cat(3,data_trials.dict); 
% rescale for each atom - to be deleted
for iAtom = 1:nAtoms
    ID(:,:,iAtom) = rescale(ID(:,:,iAtom),0,1);
end
% ID = reshape(ID,nChannels,nTimes,nAtoms);
ID = permute(ID,[1,3,2]);
nCols = ceil(nAtoms/nRows);
tcl=tiledlayout(1,nCols, 'Padding', 'tight', 'TileSpacing', 'tight');

for iCol = 1:nCols
    nexttile(tcl);
    row_start   = (iCol-1) * nRows+1;
    row_end     = iCol * nRows;
    idxs        = row_start:row_end;
    idxalt      = row_start:nAtoms; % last row
    if length(idxalt)<length(idxs)
        idxs    = idxalt;
    end
    ImgID       = ID(:,idxs,:);
    nPlots      = length(idxs);
    ImgID       = reshape(ImgID,nChannels*nPlots,nTimes);

    imagesc(1-ImgID); colormap gray; axis on; set(gca,'YDir','reverse');
    % plot names
    for il=1:nPlots
        yline(ylinesup(il)-1,'--',labnames{idxs(il)});
    end
    % plot colors
    istart  = 0;
    hold on;
    % t_start  = 1;
    t_start  = nTimes;
    t_end    = nTimes+(nTimes/4);
    for iy = 1:nPlots
        iend = ylinesup(iy);
        H(iy)=fill([t_start,t_start,t_end,t_end],[istart,iend,iend,istart],cmaps(idxs(iy),:));
        istart=iend;
        set(H(iy),'facealpha',0.08)
    end
    yline(nTimes);
    % xlim([t_start,t_end])
    xlim([times(1),t_end])
    yticks([]);
    xlabel('time')
    xt  =xticks;
    xticks(xt(xt<times(end)))
end

if nargout>0
    h=hfig;
end
%% execinfo
if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end