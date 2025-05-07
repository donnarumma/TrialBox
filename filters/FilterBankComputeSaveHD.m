function   [data, out] = FilterBankComputeSaveHD(data,par)
% function [data, out] = FilterBankComputeSaveHD(data,par)

execinfo=par.exec;
if ~isempty(execinfo); fprintf('Function: %s ',mfilename);t=tic; end
xfld        = 'time';
InField     = par.InField;
OutField    = par.OutField;
for in=1:length(data)
    if strcmp(par.FilterBank,'One')
        par.hd1 = par.hd(in).hd;
        [data(in).(OutField),hd] = FiltBankOneSaveHD(data(in).(InField), par);
    elseif strcmp(par.FilterBank,'EEGbands')
        data(in).(OutField) = FiltBankEEGbands(data(in).(InField), par);
    elseif strcmp(par.FilterBank,'Nine')
        data(in).(OutField) = FiltBankNine(data(in).(InField), par);
    elseif strcmp(par.FilterBank,'Prior')
        data(in).(OutField) = FiltBankPrior(data(in).(InField), par);
    end
    data(in).([xfld OutField])=data(in).([xfld InField]);
    data(in).hd = hd;
end

if ~isempty(execinfo); out.exectime=toc(t); fprintf('| Time Elapsed: %.2f s\n',out.exectime); end