function params = getDefaultLFPparams()
% function params = getDefaultLFPparams()

params.RS_Factor        =20; % resampling factor    
% apply filter params
params.filterOrder      =100;
% NOTCH params
params.NOTCHfilter      =0;  % apply NOTCH filter
params.NOTCHfilterLow   =45; % low   NOTCH freq
params.NOTCHfilterHigh  =55; % high  NOTCH freq

params.LPfilter         =0;  % apply   LP  filter
params.LPfilterLow      =1;  % low     LP  freq
params.LPfilterHigh     =100;% high    LP  freq

params.SubtractMean     =0;  % subtract mean values

params.SubtractDCoffset =1;  % subtract DC OFFSET
% remove EYE artifacts
params.removeEYE        =0;

params.removeRAW        =1;   % remove field RAW after extraction