function   V = cspProjection(eeg, W)
% function V = cspProjection(eeg, W)

% This function performs the decoding of spatial filtering using the CSP algorithm.
%   INPUT:
%   'eeg' is an array with EEG filtered at differenumTrial pass-bands;
%       dimensions are numb_channel•trials*samples(TW)*number_of_filters
%   'ch' are the vector channels to consider for the eeg signal extraction;
%   W{ic,1} are the projection matrices associated to the ic class
%       dimensions are numb_channel x num_of_consider_features
%       the num_of_consider_features are at least numb_channel if m≠0 is
%       2*m
%   OUTPUT:
%   V{i,1} are the features extracted from 'eeg' through W_all{i,1} projection
%   matrices
%
%   UPDATE: 2024/01/29

W = cellfun(@(x) cat(3, x{:}), num2cell(W, 2), 'UniformOutput', false);

numFilter   = size(W{1,1},3);
numClass    = size(W,1);
numTrial    = size(eeg,4);
V           = cell(numClass,1);

eeg = permute(eeg,[1 2 4 3]);
numcomponents=size(W{1},2)/2;
% Features
for i = 1:numFilter
    for j = 1:numTrial
        E = eeg(:,:,j,i);
        for k = 1:numClass
            C_class = (W{k,1}(:,:,i)'*E)*E'*W{k,1}(:,:,i); %% warning parentesi
            V{k,1}(j,(i-1)*2*numcomponents+1:i*2*numcomponents) = transpose(log(diag(C_class)/trace(C_class)));
        end
    end
end
end