function   [dataTrials,out] =pcaSynergyEncode(dataTrials,par)
% function dataTrials =pcaSynergyEncode(dataTrials,par)
execinfo=par.exec;
if ~isempty(execinfo); t=tic; end

InField       = par.InField;
OutField      = par.OutField;

xfld = 'time';

if length(size(dataTrials(1).(InField))) > 2
    data_3d = cat(4,dataTrials.(InField));
    data_3d = permute(data_3d,[4,1,2,3]);
    data_3d  = reshape(data_3d,size(data_3d,1),size(data_3d,2)*size(data_3d,3)*size(data_3d,4)); % trials x variables*time
else
    data_3d   = cat(3,dataTrials.(InField));                                      % variables x time x trials.
    data_3d   = permute(data_3d,[3,1,2]);                                                         % trials x variables x time
    data_3d   = reshape(data_3d,size(data_3d,1),size(data_3d,2)*size(data_3d,3));     % trials x variables*time
end
Zpca_test=data_3d*par.Wpca;
for iT = 1:size(data_3d,1)
    dataTrials(iT).(par.OutField) = Zpca_test(iT,:);
    dataTrials(iT).([xfld OutField])=dataTrials(iT).([xfld InField])(end);
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('Function: %s | Time Elapsed: %.2f s\n',mfilename,out.exectime); end