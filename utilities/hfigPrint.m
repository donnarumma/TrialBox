function   hfigPrint(hfg,par)
% function hfigPrint(hfg,par)

nfields             = fieldnames(hfg);
OPTIONS             = par.OPTIONS;
save_dir            = par.save_dir;
sub_dir             = par.sub_dir;
pdf_file            = par.pdf_file;
SEP                 = filesep;
PDF_dir             = [];   

for ifld = 1:length(nfields)
    fname  = nfields{ifld};
    hfig_e = hfg.(fname);
    for iopt=1:length(OPTIONS)
        t   hfigPrint= tic;
        if sub_dir
            add=upper(OPTIONS{iopt});
        else
            add='';
        end
        IMG_dir=[save_dir add SEP];
        if ~isfolder(IMG_dir);    mkdir(IMG_dir);    end
        opt = lower(OPTIONS{iopt});
        sv  = [IMG_dir fname '.' opt];
        if ismember(opt,{'png'})
            exportgraphics(hfig_e,sv);
        elseif ismember(opt,{'fig'})
            savefig(hfig_e,sv);
        else
            exportgraphics(hfig_e,sv,'ContentType','vector');
        end
        if strcmp(opt,'pdf')
            PDF_dir=IMG_dir;
        end
        fprintf('Saved %s. Elapsed time %gs\n', sv,toc(t));
    end
end
if  ~isempty(PDF_dir) & ~isempty(pdf_file)
    %% compile latex in a temporary directory
    texparams           = latexCompile_ImageDirectoryParams();
    texparams.wm        = par.wm;
    if ~isempty(par.ltx_dir)
        texparams.ltx_dir   = par.ltx_dir;
    end
    [~,outnf]   = latexCompile_ImageDirectory(PDF_dir,texparams);
    %% rename and move compiled pdf file
    fprintf('Compiled %s\n',pdf_file)
    movefile(outnf,pdf_file);
end
return
