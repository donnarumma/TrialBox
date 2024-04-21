function   par = latexCompile_ImageDirectoryParams(parSource)
% function par = latexCompile_ImageDirectoryParams(parSource)
SEP             = filesep;
par.exec        = true;
par.wm          = 200;
par.nSubfigures = 1;
par.file_type   ='pdf';
par.ltx_dir     = [SEP 'tmp' SEP num2str(randi(10^6,1,1)) SEP];
par.pdf_dir     = ['.' SEP];
try
    fnames=fieldnames(parSource);
    for in=1:length(fnames)
        par.(fnames{in})=parSource.(fnames{in});
    end
catch
end
end