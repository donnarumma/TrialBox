function [Ind_sel,position_Ind,ind_com] = MutualInformation(data,label,sel_class,m,k)

% function Ind_final = MIBIF_decode(V,label,m,k)
% MIBIV: Mutual Information-Based Best Individual Feature (MIBIV)
% This function evaluates the index of the best features for each class.
% INPUT:
% 'data': a matrix containing on the row the trials and on the column the
% features
% 'label': a vector containing the class labels
% 'm': the number elements to consider
% 'k': the maximum number of mutual information elements to include as features
% OUTPUT:
% 'Ind_final': the final index containing the indices of the projection matrix features to consider for classification.


numTrial    = size(label,1);
len_ci = sum(label == sel_class);
len_Vi = size(data,2);

% enumTrialropy of class vs not class
P = len_ci/numTrial;
H = -(P*log2(P) + (1-P)*log2(1-P));
I = zeros(len_Vi,1);
for j = 1:len_Vi
    Hc = 0;
    Hn = 0;
    for i = 1:numTrial
        % Parzen Window
        % sigma of the distribution
        sig = std(data(i,j) - data(:,j));
        h = sig*(4/(3*numTrial))^(1/5);
        % Phi
        phi = exp(-((data(i,j) - data(:,j)).^2)/(2*h^2))/sqrt(2*pi);

        % estimation of conditional probability with Parzen Window
        cp = sum(phi(label == sel_class))/len_ci;
        np = sum(phi(label ~= sel_class))/(numTrial - len_ci);

        % probability of the feature (j,i)
        p = cp*P + np*(1-P);

        % conditional probabilities
        pc = cp*P/p;
        pn = np*(1-P)/p;

        % conditional enumTrialropy for 'class' and 'not class'
        Hc = Hc + pc*log2(pc);
        Hn = Hn + pn*log2(pn);
    end

    % Conditional enumTrialropy
    Hcond = -(Hc + Hn);

    % mutual information of features for desired class
    I(j) = H - Hcond; % information content (2*m*numerodifiltridelVB,1)
end

% Features selection
[~, ind] = sort(I); %ind dimensions: num_filter*num_channel

% selected indexes (first k features in descending order)
ind_sel = ind(end-(k-1):end); %Starting from the bottom, it takes the K_MIBIV features with the highest informative content.
% complementary indexes (according to CSP)
if (m == 0)
    ind_com = ind_sel;
else
    ind_com = (4*m)*ceil(ind_sel/(2*m))-(2*m-1)-ind_sel;
end

% final indexes
[Ind_sel,position_Ind] = unique([ind_sel ind_com]);
end
