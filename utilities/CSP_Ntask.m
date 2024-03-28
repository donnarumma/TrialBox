function W = CSP_Ntask(X,m)
% CSP creates the projection matrices "Wi", associated to the 4 classes,
% starting from the spectrogram matrices "Xi" separated per class.
%
%   INPUT:
%   'X' is the cell containing the 3D array compose by the signals of the
%   all classes
%
%   'm' is the number of first and last rows to consider in the projection
%   matrices 'Wi'.
%
%   OUTPUT:
%   'W' is the cell containing the projection matrix associated to the all classes;
%
%
%  authors:         M. Frosolone
%  correspondence:  m.frosolone@istc.cnr.it
%  last update:     2023/10/13


    % Covariance matrices
    numClass        = size(X,1); % number of classes
    nch             = size(X{1},1); % number of channels
    C_final         = cell(numClass,1);
    
    for cl = 1:numClass
        % Covariance of i class
        el_class        = size(X{cl},3);
        C_class         = zeros(nch,nch);
        for i = 1:el_class
            M           = X{cl}(:,:,i)*X{cl}(:,:,i)';
            C_class     = C_class + M/trace(M);
        end
        C_class         = C_class/el_class;
        C_final{cl,1}   = C_class;
    end
    
    % Composite Covariance Matrix
    C_compos            = zeros(nch,nch);
    for i = 1:length(C_final)
        C_compos        = C_compos + C_final{i,1};
    end
    
    % Projection matrices
    W               = cell(numClass,1);
    A               = cell(numClass,1);
    for cl = 1:numClass
        % Complete projection matrices
        [W_class, A_class] = eig(C_final{cl,1},C_compos);
        W{cl,1}            = W_class;
        A{cl,1}            = A_class;
    end
    
    % Sorting
    for s = 1:numClass
        [~, ind_class] = sort(diag(A{s,1}));
        W{s,1}         = W{s,1}(:,ind_class);
    end
    
    if m ~= 0
        for i = 1:numClass
            W{i,1}     = [W{i,1}(:,1:m), W{i,1}(:,end-m+1:end)];
        end
    end

end

