function   [W_Filter, out] = CSPDictionary(X_in, label_classes, cspcomponents,par)
% function   [W_Filter, out] = CSPDictionary(X_in, label_classes, cspcomponents,par)
% This function get the projection matrix WCSP following the CSP method.
%
%   INPUT:
%   'X_in' is a cell with Data. in each cell X_in{it} space dimension (e.g.
%   number of channel) per time dimension per number of filters of FILTER
%   BANK - FB
%   'cspcomponents' is the number of CSP components to consider.
%
%   OUTPUT:
%   'W_Filter' is a cell containing the CSP projection matrix related to all classes and all filters.
%   for each class ic the relative projection matrix is WCSP{ic};
%
%   UPDATE: 2024/01/19
%
execinfo=par.exec;
if ~isempty(execinfo); t=tic; end

typeClass       = unique(label_classes);
numClass        = length(typeClass);
numFilter       = size(X_in,3);
X_class         = cell(numClass,numFilter);

for ifil=1:numFilter
    X_fil = X_in(:,:,ifil,:);
    for ic=1:numClass
        % X_class{ic,1}   = cat(3,X_in{label_classes==typeClass(ic)});
        % X_class{ic,1}   = cat(3,X_in{label_classes==typeClass(ic)});
        X_class{ic,ifil}   = X_fil(:,:,label_classes==typeClass(ic));
    end
end

W_Filter               = cell(numClass,numFilter);
% projection matrices
for idf = 1:numFilter
    X_cf                    = X_class(:,idf);
    W_Filter(:,idf)    = CSP_Ntask(X_cf,cspcomponents);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end
out.W      = W_Filter;