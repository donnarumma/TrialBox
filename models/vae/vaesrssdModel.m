function   [data_trials, out]=vaesrssdModel(data_trials,par)
% function [data_trials, out]=vaesrssdModel(data_trials,par)
execinfo                = par.exec;
if ~isempty(execinfo); t=tic; fprintf('Function: %s ',mfilename); end

InField                 = par.InField;
% cat(3,dataTrials.(fn)) put trials on the third dimension,
X3d                     = cat(3,data_trials.(InField));                     % nCells x nTimes x nTrials
% X3d         = permute(X3d,[3,1,2]);                             % nTrials x nCells x nTimes
% each row is a "synergy" -> [x_1(1), ... x_nChannels(1), x_1(2), ... x_nChannels(2), ... , x_1(nTimes), ... x_nChannels(nTimes),    

X4d                     = reshape(X3d,size(X3d,1),size(X3d,2),1,size(X3d,3)); % nTrials x nChannels x 1 x nTimes
[nCells, nTimes, nChannels, nTrials] = size(X4d);
if canUseGPU
    fprintf('using GPU ');
    X4d = gpuArray(X4d);
end
nValids                 = round(nTrials*par.validPercent/100);
iTrials                 = randperm(nTrials);
iValid                  = iTrials(1:nValids);
iTrain                  = iTrials(nValids+1:end);
X4dTrain                = X4d(:,:,:,iTrain);
X4dValid                = X4d(:,:,:,iValid);
if isempty(X4dValid)
    validationFrequency = false;
else
    validationFrequency = par.validationFrequency;
    nTrials             = length(iTrain);
    dsValid             = arrayDatastore(X4dValid,IterationDimension=4);
    mbqValid            = minibatchqueue(                                           ...
                                        dsValid,                                    ...
                                        par.MBnumOutputs,                           ...
                                        MiniBatchSize       = par.miniBatchSize,    ...
                                        MiniBatchFcn        = @preprocessMiniBatch, ...
                                        MiniBatchFormat     = "SSCB",               ...
                                        PartialMiniBatch    = "discard");

end

dsTrain                 = arrayDatastore(X4dTrain,IterationDimension=4);
mbqTrain                = minibatchqueue(                                           ...
                                        dsTrain,                                    ...
                                        par.MBnumOutputs,                           ...
                                        MiniBatchSize       = par.miniBatchSize,    ...
                                        MiniBatchFcn        = @preprocessMiniBatch, ...
                                        MiniBatchFormat     = "SSCB",               ...
                                        PartialMiniBatch    = "discard");


% encoder net
if isempty(par.netE)
    imageSize               = [nCells, nTimes, nChannels];
    layersE                 = [ 
                                imageInputLayer(imageSize,Normalization="none") ,                                       ...
                                convolution2dLayer(par.kernsize,  par.numLatentChannels,Padding="same",Stride=2),       ...
                                reluLayer,                                                                              ...                   
                                convolution2dLayer(par.kernsize,2*par.numLatentChannels,Padding="same",Stride=2),       ...
                                reluLayer,                                                                              ...
                                fullyConnectedLayer(2*par.numLatentChannels),                                           ...
                                samplingLayer];
    netE                    = dlnetwork(layersE);
else
    netE                    = par.netE;
end
                        
% decoder net
if isempty(par.netD)
    % projectionSize          = [par.projsize, par.projsize, 2*par.numLatentChannels];
    % layersD                 = [
    %                             featureInputLayer(par.numLatentChannels),                                               ...
    %                             projectAndReshapeLayer(projectionSize),                                                 ...
    %                             transposedConv2dLayer(par.kernsize,2*par.numLatentChannels,Cropping="same",Stride=2),   ...
    %                             reluLayer,                                                                              ...
    %                             transposedConv2dLayer(par.kernsize,  par.numLatentChannels,Cropping="same",Stride=2),   ...
    %                             reluLayer,                                                                              ...
    %                             transposedConv2dLayer(par.kernsize,nChannels,Cropping="same"),                          ...
    %                             sigmoidLayer];
    % projectionSize          = [par.projsize, par.projsize, 2*par.numLatentChannels];
    layersD                 = [
                                featureInputLayer(par.numLatentChannels), ...
                                projectionLayer(imageSize),              ...projectAndReshapeLayer(imageSize), ... 
                                ];
    netD                    = dlnetwork(layersD);
else
    netD                    = par.netD;
end

numIterationsPerEpoch   = ceil(nTrials / par.miniBatchSize);
numIterations           = par.numEpochs * numIterationsPerEpoch;

if par.graphplot
    monitor             = trainingProgressMonitor(  ...
                                                    Metrics             = "Loss",   ...
                                                    Info                = "Epoch",  ...
                                                    XLabel              = "Iteration");
else
    monitor.Stop        = false;
end
% learning
trailingAvgE            = [];
trailingAvgSqE          = [];
trailingAvgD            = [];
trailingAvgSqD          = [];
epoch                   = 0;
iteration               = 0;

% Loop over epochs.
while epoch < par.numEpochs && ~monitor.Stop
    epoch = epoch + 1;
    
    % Shuffle data.
    shuffle(mbqTrain);
    
    % Loop over mini-batches.
    while hasdata(mbqTrain) && ~monitor.Stop
        iteration = iteration + 1;

        % Read mini-batch of data
        Xbatch = next(mbqTrain);

        % Evaluate loss and gradients
        [loss,gradientsE,gradientsD] = dlfeval(@modelLoss,netE,netD,Xbatch,par);
        
        % Update learnable parameters
        [netE,trailingAvgE,trailingAvgSqE]      = adamupdate(netE, ...
            gradientsE,trailingAvgE,trailingAvgSqE,iteration,par.learnRate);
        [netD, trailingAvgD, trailingAvgSqD]    = adamupdate(netD, ...
            gradientsD,trailingAvgD,trailingAvgSqD,iteration,par.learnRate);

        % Every validationFrequency iterations, display batch of generated images using the
        % held-out generator input
        if validationFrequency && (mod(iteration,validationFrequency) == 0 || iteration == 1) 

            % evaluate on the validation
            XValidbatch     = next(mbqValid);
            validloss       = dlfeval(@modelLoss,netE,netD,XValidbatch);
        end
        
        if par.graphplot
            % Update the training progress monitor. 
            recordMetrics(monitor,iteration,Loss=loss);
            updateInfo(monitor,Epoch=epoch + " of " + par.numEpochs);
            
            monitor.Progress    = 100*iteration/numIterations;
        end
    end
    if ~par.graphplot
        fprintf('Epoch %g: it=%g/%g | Loss: %g\n',epoch,iteration,numIterations,loss);
    end
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end
out.netE    = netE;
out.netD    = netD;
end

function   X = preprocessMiniBatch(dataX)
% function X = preprocessMiniBatch(dataX)
    % Concatenate.
    X                       = cat(4,dataX{:});
end

function   [loss,gradientsE,gradientsD] = modelLoss(netE,netD,Xdata,par)
% function [loss,gradientsE,gradientsD] = modelLoss(netE,netD,Xdata)
    % Forward through encoder.
    [Zdata,mu,logSigmaSq]   = forward(netE,Xdata);
    
    % Forward through decoder.
    Ypred                   = forward(netD,Zdata);
    
    % Calculate loss and gradients.
    loss                    = elboLossSrssd(Ypred,Xdata,Zdata,mu,logSigmaSq, ...
                                            netD.Layers(2).Weights,par);
    [gradientsE,gradientsD] = dlgradient(loss,netE.Learnables,netD.Learnables);
end

function   loss = elboLossSrssd( Ypred, Xbatch, Zdata, mu, logSigmaSq, ...
                                 Dsrssd, par)
% function loss = elboLoss(Ypred,Xbatch,mu,logSigmaSq)
    % reconstruction in the original space
    % nAtoms % Number of dictionary elements
    [nCells, nTimes,~, nTrials]  = size(Ypred);
    spG_squared         = par.spG.^2;
    % Xbatch_e = squeeze(extractdata(Xbatch)); Ypred_e = squeeze(extractdata(Ypred)); Zdata_e             = extractdata(Zdata); 
    % logSigmaSq_e = extractdata(logSigmaSq); mu_e = extractdata(mu);
    
    % loss1       = norm( Xbatch_e - Ypred_e, 'fro' )^2 / norm( Xbatch_e, 'fro' )^2;
    loss1       = mse(Ypred,Xbatch);
    % KL divergence.
    KL          = -0.5 * sum(1 + logSigmaSq - mu.^2 - exp(logSigmaSq),1);
    KL          = mean(KL);
    loss1       = loss1 + KL;
    
    % reconstruction in the manifold space
    % loss2      = norm( U - X*C', 'fro' ) ^2 / (n*r)*par.eta;
    
    % Omega_approx
    [inv_Zeta, Eta] = mexUpdateEta( reshape(Dsrssd,[nCells*nTimes,par.nAtoms]), spG_squared, par );
    norm_eta        = sum( sum( Eta .^ 2, 1, 'omitnan' ) .^ (1/2) ); 
    quadratic_term  = sum( sum( reshape(Dsrssd.^2, [nCells*nTimes,par.nAtoms]) .* inv_Zeta, 1,'omitnan') );   
    loss3           = 0.5 * quadratic_term + 0.5 * norm_eta;
    
    % loss            = loss1 + par.lambda *loss3;
    % return
   
    % constraint on manifold
    loss4           = sum(abs(Zdata(:)))/(nTrials*par.nAtoms); % 'OmegaU'
     
    % Combined loss
    loss            = loss1;                    % ? + classical value KL
    % loss            = loss + loss2;
    loss            = loss + par.lambda *loss3; % structure
    loss            = loss + par.lambda2*loss4; % ? sparsity on Zdata? necessity of abs?
    fprintf('Loss1: %g | Loss3: %g | Loss4 %g | Loss3*lamba: %g | Loss4*lamba: %g\n',loss1, loss3, loss4,loss3*par.lambda,loss4*par.lambda2);
    return
    % loss            = dlarray(loss,'CB');
end