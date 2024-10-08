function   W_Filter = cspDictionary(X_in, label_classes, cspcomponents)
% function W_Filter = cspDictionary(X_in, label_classes, cspcomponents)
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
typeClass       = unique(label_classes);
numClass        = length(typeClass);
numFilter       = size(X_in,3);
X_class         = cell(numClass,numFilter);

for ifil=1:numFilter
    X_fil = X_in(:,:,ifil,:);
    for ic=1:numClass
        X_class{ic,ifil}   = X_fil(:,:,label_classes==typeClass(ic));
    end
end

W_Filter               = cell(numClass,numFilter);
% projection matrices
for idf = 1:numFilter
    X_cf               = X_class(:,idf);
    W_Filter(:,idf)    = cspNtask(X_cf,cspcomponents);
end