function   dict_trials = arrangeDictionary (dict_trials,par)
% function dict_trials = arrangeDictionary (dict_trials,par)
if isempty(dict_trials)
    dict_trials=struct;
end
dictionary= par.dictionary;
nChannels = par.nChannels;
nTimes    = par.nTimes;
OutField  = par.OutField;
xfld      = par.xfld;
if isempty(par.times)
    times=1:nTimes;
else
    times=par.times;
end
[nAtoms,nChannelsxnTimes]=size(dictionary);
if nChannelsxnTimes ~= nChannels*nTimes
    fprintf('Dimensions are not consistent: in dictionary: %g, in params: %g x %g = %g\n',nChannelsxnTimes,nChannels,nTimes,nChannels*nTimes)
    return;
end
for iAtom=1:nAtoms
    dict_trials(iAtom).(OutField)=reshape(dictionary(iAtom,:),nChannels,nTimes);
    dict_trials(iAtom).([xfld OutField])=times;
end
