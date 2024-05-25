function   Model=pcaAtom(Model,index)
% function Model=pcaAtom(Model,index)
if nargin  < 2
    Model = size(Model.Wpca,2); % return number of components if no index is given
else
    Model.Wpca=Model.Wpca(:,index); % Dictionary is nTimes x nComponents
end
