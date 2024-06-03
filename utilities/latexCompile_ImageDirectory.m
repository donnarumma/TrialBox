function   [s,pdf_out]=latexCompile_ImageDirectory(source_dir,params)
% function [s,pdf_out]=latexCompile_ImageDirectory(source_dir,params)
SEP             = filesep;
try
    wm          = params.wm;
    nSubfigures = params.nSubfigures;
    file_type   = params.file_type;
    ltx_dir     = params.ltx_dir;
catch
    wm          = 100;
    nSubfigures = 4;
    file_type   ='pdf';
    ltx_dir     = [SEP 'tmp' SEP num2str(randi(10^6,1,1)) SEP];
end
imgfiles        = dir([source_dir SEP '*.' file_type]);
nFiles          = length(imgfiles);

s=sprintf('\\documentclass[english]{article}\n\\usepackage[T1]{fontenc}\n\\usepackage[latin9]{inputenc}\n\\usepackage{caption}\n\\usepackage{float}\n\\usepackage{graphicx}\n\\usepackage{babel}\n\\usepackage{subcaption}\n\\extrafloats{200}\n\\maxdeadcycles=200\n\\usepackage{epstopdf}\n\\usepackage{comment}\n\\begin{document}\n');

%% loop 

nFigures=ceil(nFiles/nSubfigures);

for n_f=1:nFigures
    group_ind=(1:nSubfigures)+(n_f-1)*nSubfigures;
    s=sprintf('%s\n\\begin{figure}\n\\begin{centering}\\vspace{-3em}\n',s);
    for i=group_ind
        if  mod(i-group_ind(1),2)
            hg=-1;
        else
            hg=-10;
        end
        s=sprintf('%s\\hspace*{%gem}\\subcaptionbox[\\label{Fig:N%g}]{\\includegraphics[width=%gmm]{%s}}%%\n',s,hg,i,wm,[regexprep(imgfiles(i).folder,'\','/') '/' imgfiles(i).name(1:end-4)]);
        if  ~mod(i-group_ind(1),2)    
        else
            s=sprintf('%s\n',s);
        end
        if i==nFiles
            break;
        end
    end
    s=sprintf('%s\\par\\end{centering}\n\\caption{}\\label{Fig:%g}\n',s,n_f);
    s=sprintf('%s\\end{figure}\n\n',s);
end

%% close
s=sprintf('%s\\end{document}\n',s);

%% save tex file
% susbstitute unallowed chars with '_'
unallowed = {SEP '/' '\' '~'};
filename = source_dir;
for iu=1:length(unallowed)
    filename = strrep(filename,unallowed{iu},'_');
end
% avoid start and end with '_'
if strcmp(filename(1),'_')
    filename=filename(2:end);
end
if strcmp(filename(end),'_')
    filename=filename(1:end-1);
end
ltxfilename=sprintf('%s.tex',filename);
% save latex file in a temporary directory
if ~isfolder(ltx_dir)
    mkdir(ltx_dir);
end
fullpathltxfilename = sprintf('%s%s',ltx_dir,ltxfilename);
ltxfileid           = fopen(fullpathltxfilename,'w');
fprintf('Compiling %s\n',fullpathltxfilename)
fprintf(ltxfileid,'%s',s);
fclose(ltxfileid);

%% compile to pdf
curdir   = cd;
cd (ltx_dir);
% latex command -> system(latexmk -pdf -shell-escape -synctex=1 namefile.tex');
comm     = ['latexmk -pdf -shell-escape -synctex=1 ' ltxfilename];
fprintf('executing %s\n',comm);
system(comm);
cd(curdir);
% compiled pdf in pdf_out
pdf_out  = [ltx_dir filename '.pdf'];
fprintf('Created %s.pdf\n',[ltx_dir filename]);

% copy compiled pdf in out dir
% fprintf('File ready in %s.pdf\n',[pdf_dir SEP filename]);
% outnf   = [pdf_dir SEP filename '.pdf'];
% copyfile([fullpathltxfilename(1:end-4) '.pdf'], pdf_dir);