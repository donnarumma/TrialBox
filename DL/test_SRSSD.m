%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% load MNIST in Matlab
[XTrain,YTrain,anglesTrain] = digitTrain4DArrayData;
[XTest, YTest, anglesTest ] = digitTest4DArrayData;
X3dTrain                    = squeeze(XTrain);
X3dTest                     = squeeze(XTest);

%% arrange
[nChannels,nTimes,nTrials]  = size(X3dTrain);
% nChannels*nTimes x nTrials 
XdataTrain                  = reshape(X3dTrain,nChannels*nTimes,nTrials);
XdataTest                   = reshape(X3dTest,nChannels*nTimes,nTrials);
% nTrials x nChannels*nTimes
XdataTrain                  = XdataTrain';
XdataTest                   = XdataTest';

%% Test settings srssd
%--------------------------------------------------------------------------
logLambdaA          =-7;%[-8,-6];
% lambdaAtom = 10^-5;
% lambdaCoeff = 0.5;
logLambdaA          =-5;%[-8,-6];

lambdaARange        = 10.^logLambdaA;
lambdaCRange        = 1.5;%0.9;%0.5; % 0.02,0.8
nAtoms              = 100;
selectedStruct      = 'ContinuousTemporal';
runFields           = {'lambda','lambda2'};
% spG                 = setActionSPG(nTimes,nChannels,selectedStruct);
mode.rectangle      = true; 
mode.pi4            = true;
spG                 = get_groups( mode, [nChannels nTimes] );

par.srssdModel      = srssdCreateLearnParams(nAtoms,     ...
                                             lambdaARange, ...
                                             lambdaCRange, ...
                                             spG,          ...
                                             runFields,    ...
                                             [],           ...
                                             []);
par.srssdModel.normparam  = 0.5;
par.srssdModel.eta        = 1; % default 0

%% model learning
% sr_ssd: n x d -> n=observations d=features -> n=nTrials, d=nChannels*nTimes
% [U, D, C]   = sr_ssd(X, spG, params);
fprintf('Learning srssd Model...'); t=tic; 
[Manifold, Dsrssd, Wsrssd]  = sr_ssd(XdataTrain, spG, par.srssdModel);
fprintf('Elapsed time %.2f s\n',toc(t));

% Dsrssd -> nChannels*nTimes x nAtoms
% Wsrssd -> nAtoms x nChannels*nTimes  ????

%% manifold construct by sim_sr_ssd    -> see below for an alternative
% nTrials x nAtoms
fprintf('Projecting ...'); t=tic;
ZdataTrain                  = sim_sr_ssd(XdataTrain, par.srssdModel, Dsrssd);
ZdataTest                   = sim_sr_ssd(XdataTest, par.srssdModel, Dsrssd);
fprintf('Elapsed time %.2f s\n',toc(t));
%% manifold construct by Coding Matrix -> see above for an alternative
% nTrials x nAtoms
ZdataTrain                  = XdataTrain*Wsrssd'; % XdataTrainrec   = ((XdataTrain*C')*D');
% nTrials x nAtoms 
ZdataTest                   = XdataTest*Wsrssd';  % XdataTestrec    = ((XdataTest*C')*D');
%% reconstruction
% nTrials x nChannels*nTimes
XdataTrainrec               = ZdataTrain*Dsrssd';
XdataTestrec                = ZdataTest*Dsrssd';
% reconstruction error
errorXTrain                 = RMSError(XdataTrainrec,XdataTrain);
errorXTest                  = RMSError(XdataTestrec,XdataTest);
out.errorXTrain             = errorXTrain;
out.errorXTest              = errorXTest;
disp(out)
%% Learn classifier
t=tic;
fprintf('Learning Classifier...')
par.kfold           = 4;
par.numIterations   = 300;
cvp                 = cvpartition(YTrain,'kfold',par.kfold,'Stratify',true);

% fitdiscr input -> nTrials x nAtoms
mdl = fitcdiscr(ZdataTrain, YTrain,...
    'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions', ...
    struct('Optimizer','randomsearch','CVPartition',cvp,'MaxObjectiveEvaluations',par.numIterations, 'AcquisitionFunctionName','expected-improvement-plus','ShowPlots',false,'Verbose',0));
fprintf('Elapsed time %g s\n',toc(t));
%% accuracy on Test
[YTest_pred,probs]= predict(mdl, ZdataTest); % nTrials x nClasses
out.accuracy      = mean(YTest_pred==YTest);
out.accuracyclass = accuracy4classes(YTest_pred,YTest);
disp(out)
%% plot reconstructions
col             = [0.5,0.5,0.7]; % col borders
hmShow          = 49;
% train
numObservations = size(XTrain,4);
idx             = randperm(numObservations,hmShow);
ImgTrain        = reshape (XdataTrain(idx,:)',nChannels,nTimes,hmShow);
ImgTrainRec     = reshape (XdataTrainrec(idx,:)',nChannels,nTimes,hmShow);
subplot(1,2,1)
montage (ImgTrain, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Original');
subplot(1,2,2);
montage (ImgTrainRec, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Reconstructed');
sgtitle(sprintf('Train (RMS: %.4f)',errorXTrain))
% test
numObservations = size(XTest,4);
idx             = randperm(numObservations,hmShow);
ImgTest         = reshape (XdataTest(idx,:)',nChannels,nTimes,hmShow);
ImgTestRec      = reshape (XdataTestrec(idx,:)',nChannels,nTimes,hmShow);
figure;
subplot(1,2,1)
montage (ImgTest, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Original');
subplot(1,2,2);
montage (ImgTestRec, 'BorderSize', [1,1], 'BackgroundColor', col);
title('Reconstructed');
sgtitle(sprintf('Test (RMS: %.4f)',errorXTest))
%% plot dictionaries 
% encoding
figure;
tiledlayout(1,2, 'Padding', 'tight', 'TileSpacing', 'tight');
nexttile;
hold on;
title ('Encoding Dictionary')
IW = Wsrssd';
% IW = Wsrssd; % consider the space of Atoms and not figures
% IW = rescale(IW,0,1); % uniformly scale
for iAtom = 1:size(IW,2)
    IW(:,iAtom) = rescale(IW(:,iAtom),0,1);
end
IW=reshape (IW,nChannels,nTimes,nAtoms);
montage (IW, 'BorderSize', [2,2], 'BackgroundColor', col);
% decoding
nexttile;
hold on;
title ('Decoding Dictionary')
ID = Dsrssd;
% ID = rescale(ID,0,1); % uniformly scale
for iAtom = 1:nAtoms
    ID(:,iAtom) = rescale(ID(:,iAtom),0,1);
end
ID=reshape (ID,nChannels,nTimes,nAtoms);
montage (ID, 'BorderSize', [2,2], 'BackgroundColor', col);
return

%% plot dictionaries 
% encoding
hfig=figure;
% nexttile;
tcl=tiledlayout(ceil(sqrt(nAtoms)),ceil(sqrt(nAtoms)), 'Padding', 'tight', 'TileSpacing', 'tight');
% tcl=tiledlayout(nAtoms,1, 'Padding', 'tight', 'TileSpacing', 'tight');
han                 = axes(hfig,'visible','off'); 
title ('Encoding Dictionary')
IW = Wsrssd';
% IW = rescale(IW,0,1); % uniformly scale
for iAtom = 1:nAtoms
    nexttile(tcl)
    IW(:,iAtom) = rescale(IW(:,iAtom),0,1);
    imagesc(reshape(IW(:,iAtom),nChannels,nTimes)); colormap gray
    % axis off;
end
%%
IW=reshape (IW,nChannels,nTimes,nAtoms);
montage (IW, 'BorderSize', [2,2], 'BackgroundColor', col);%,'Size',[nAtoms, 1]);
%%
% decoding
nexttile;
hold on;
title ('Decoding Dictionary')
ID = Dsrssd;
% ID = rescale(ID,0,1); % uniformly scale
for iAtom = 1:nAtoms
    ID(:,iAtom) = rescale(ID(:,iAtom),0,1);
end
ID=reshape (ID,nChannels,nTimes,nAtoms);
montage (ID, 'BorderSize', [2,2], 'BackgroundColor', col);
return
%%
figure;
IW = Wsrssd';
imagesc(reshape(permute(IW,[3,1,2]),nAtoms*nChannels,nTimes));colormap gray;
%%
figure;
ID = Dsrssd;
ciccio=[];
for iAtom = 1:nAtoms
    ID(:,iAtom) = rescale(ID(:,iAtom),0,1);
    ciccio = [ciccio;reshape(ID(:,iAtom),nChannels,nTimes)];
end
% ciccio = reshape(ciccio,nAtoms,nChannels,nTimes);
imagesc(ciccio);colormap gray; axis on;
%%
figure;
nAtoms      =100;
ylines      =1:nChannels:nChannels*nAtoms;
ylinesup    =ylines(2:end);
ylinedown   =ylines(1:end-1);
for iAtom=1:nAtoms
    labnames{iAtom}=sprintf('A%g',iAtom);
end
idxToPlot   = 1:20;
cmaps       =linspecer(length(idxToPlot));

ID          = Dsrssd(:,idxToPlot);
nAtoms      = length(idxToPlot);
labnames    = labnames(idxToPlot);
for iAtom = 1:nAtoms
    ID(:,iAtom) = rescale(ID(:,iAtom),0,1);
end
%
ID = reshape(ID,nChannels,nTimes,nAtoms);
ID = permute(ID,[1,3,2]);
ID = reshape(ID,nChannels*nAtoms,nTimes);
imagesc(1-ID); colormap gray; axis on
%
for il=1:nAtoms
    yline(ylinesup(il)-1,'--',labnames{il});
end

istart  = 0;
hold on;
t_start  = 1;
t_end    = nTimes+nTimes/10;
for iy = 1:nAtoms
    iend = ylinesup(iy);
    H(iy)=fill([t_start,t_start,t_end,t_end],[istart,iend,iend,istart],cmaps(iy,:));
    istart=iend;
    set(H(iy),'facealpha',0.08)
    HylinesL{iy}=['Step: ' num2str(iy)];
end
xlim([t_start,t_end])
yticks([]);
%%
ID=reshape (ID,nChannels,nTimes,nAtoms);
ID=ID(:,:);
% ID = reshape(ID,[3,1,2]);

imagesc(ID');colormap gray;
