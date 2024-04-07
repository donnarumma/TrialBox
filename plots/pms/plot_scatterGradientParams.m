function   par = plot_scatterGradientParams(parSource)
% function par = plot_scatterGradientParams(parSource)
par.exec        = true;
par.hfig        = [];
par.InField     = 'y';                  
par.InGradient  = 'gradients';
par.fine        = 10;                   % round approximation of gradient colorbar
par.p           = [0,0,1500,600];       % plot size
par.fs          = 14;                   % font size
par.lats        = [1,2,3];              % directions to be plot     
par.ms          = 10;                   % marker size
par.widthsub    = 15;                   % subplot size of the main plot
par.reverse     = [false,false,false];   % reverse directions axis
par.label       = 'time';
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end