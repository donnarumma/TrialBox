function   explained = pca_explained_variance(X_data,Wpca)
% function explained = explainedVariance(X_data,par)

% par.InField='encode';
% par.OutField='decode';
% [X_data,out]            =par.encode(X_data,par.encode);
% [X_data_rec,out]        =par.decode(X_data,par.decode);
% 

% X_data=rand(10,20); % nChannels x nTimes
% [Wpca, Zpca, ~, ~, explainedpca, mupca]= pca(X_data');
mu                      = mean(X_data,2);
%  Zpca_rec             =     Xdata' * Wpca    % TxN x NxK 
%  X_data_rec           = Zpca_data  * Wpca' + mu
nComps                  = size(Wpca,2);
nTimes                  = size(X_data,2);
murep                   = repmat(mu,1,nTimes);
% Z_data                  = (X_data - murep)' * Wpca;
% X_data_rec              = Wpca * Z_data' + murep;  

explained               = nan(nComps,1);
   
    %
for iComp=1:nComps
    % X_data            -> nChannels x nTimes
    % Wpca              -> nChannels x nComps
    %
    % Z_data            -> nTimes x nComps
    % X_data' * Wpca;   -> nTimes x nChannels * nChannels x nComps
    % Z_data            = (X_data - murep)' * Wpca(:,iComp);   

    % X_data_rec        -> nChannels x nTimes
    % Wpca * Z_data'    -> (nChannels x nComps * nComps x nTimes)
    % X_data_rec        = Wpca(:,iComp) * Z_data' + murep;  

    % Z_data * Wpca'    -> (nTimes x nComps * nComps x nChannels)'
    % X_data_rec          = (Z_data * Wpca(:,iComp)')' + murep;  
 
    X_mu                = X_data-murep;
    % X_data_rec          = (X_mu' * Wpca(:,iComp) * Wpca(:,iComp)')' + murep;
    % explained(iComp)    = (1-(var(X_data(:)-X_data_rec(:)))/var(X_mu(:)))*100;
   
    X_mu_rec            = (X_mu' * Wpca(:,iComp) * Wpca(:,iComp)')'; 
    explained(iComp)    = (1-(var(X_mu(:)-X_mu_rec(:)))/var(X_mu(:)))*100;
   
    % explained(iComp)    = (1-(var(X_data(:)-X_data_rec(:)))/var(X_data(:)))*100;
   
end