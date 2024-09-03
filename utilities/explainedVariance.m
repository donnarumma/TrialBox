function   explained = explainedVariance(data_trials,par)
% function explained = explainedVariance(data_trials,par)
InField                 = par.InField;         
OutField                = 'rec';
encfield                = 'pca';

nAtoms                  = par.getAtom(par.Model);
explained               = nan(nAtoms,1);

X_data                  = [data_trials.(InField)];
par.Model.exec          = 0;
par.Model.exec          = 0;
par.Model.InField       = InField;
    
for iAtom = 1:nAtoms
    Model               = par.getAtom(par.Model,iAtom);
    
    Model.InField       = InField;
    Model.OutField      = encfield;
    data_trials_enc     = par.Encoder(data_trials,Model);

    Model.InField       = encfield;
    Model.OutField      = OutField;
    data_trials_rec     = par.Decoder(data_trials_enc,Model);

    X_data_rec          = [data_trials_rec.(OutField)];
    explained(iAtom)    = (1-(var(X_data(:)-X_data_rec(:)))/var(X_data(:)))*100;   
end
               