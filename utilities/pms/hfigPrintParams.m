function   par = hfigPrintParams(parSource)
% function par = hfigPrintParams(parSource)
SEP             = filesep;
par.exec        = true;
par.save_dir    = ['.' SEP];
par.sub_dir     = true; 
par.OPTIONS     = {'PDF','FIG'}; % e.g. par.OPTIONS = {'EPS','PNG','PDF','FIG'};
par.pdf_file    = [];
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end