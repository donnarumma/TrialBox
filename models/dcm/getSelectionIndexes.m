function  isel=getSelectionIndexes (idirs,iconds)
%function isel=getSelectionIndexes (idirs,iconds)
% idirs               = 1:8;%17;
idirs               = idirs'; 
isel                = idirs+[0*8,1*8,2*8];
isel                = isel(:,iconds);
isel                = isel(:);