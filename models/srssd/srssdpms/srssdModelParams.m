function   par = srssdModelParams(parSource)
% function par = srssdModelParams(parSource)
par.exec            = true;
par.InField         = 'y';
par.nAtoms          = 100;
par.spG             = [];
par.normparam       = 0.5;
par.r               = 100;      % dictSize -> number of Atoms
par.it0             = 10;       % Cost function displayed every 25 iterations
par.min_delta_cost  = 0.01;     % Stopping criterion: relative decrease in the cost function smaller than 0.1
par.posV            = 0;
par.posU            = 0;
par.lambda          = 10^-7;    %lambdaARange -> Regolarization with respect to atoms
par.lambda2         = 0.5;      %lambdaCRange -> Regolarization with respect to Coefficients
par.max_it          = 100;
par.max_it_U        = 10;
par.max_it_V        = 10;
par.eta             = 0;
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end