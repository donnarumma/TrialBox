function params=getLFPparamsDonnarumma(reference)

params=getDefaultLFPparams;

switch reference
    case 1
        params.NOTCHfilter      =1;  % apply NOTCH filter
        params.LPfilter         =1;  % apply   LP  filter
        params.SubtractMean     =1;  % subtract mean values
        params.SubtractDCoffset =1;  % subtract DC OFFSET
        params.removeEYE        =1;  % remove EYE artifacts
    case 2
        params.NOTCHfilter      =0;  % apply NOTCH filter
        params.LPfilter         =1;  % apply   LP  filter
        params.SubtractMean     =1;  % subtract mean values
        params.SubtractDCoffset =1;  % subtract DC OFFSET
        params.removeEYE        =1;  % remove EYE artifacts
    case 3
        params.NOTCHfilter      =1;  % apply NOTCH filter
        params.LPfilter         =0;  % apply   LP  filter
        params.SubtractMean     =1;  % subtract mean values
        params.SubtractDCoffset =1;  % subtract DC OFFSET
        params.removeEYE        =1;  % remove EYE artifacts
    case 4
        params.NOTCHfilter      =0;  % apply NOTCH filter
        params.LPfilter         =0;  % apply   LP  filter
        params.SubtractMean     =1;  % subtract mean values
        params.SubtractDCoffset =1;  % subtract DC OFFSET
        params.removeEYE        =1;  % remove EYE artifacts
    case 5
        params.NOTCHfilter      =0;  % apply NOTCH filter
        params.LPfilter         =0;  % apply   LP  filter
        params.SubtractMean     =0;  % subtract mean values
        params.SubtractDCoffset =1;  % subtract DC OFFSET
        params.removeEYE        =1;  % remove EYE artifacts
    case 6
        params.NOTCHfilter      =0;  % apply NOTCH filter
        params.LPfilter         =0;  % apply   LP  filter
        params.SubtractMean     =0;  % subtract mean values
        params.SubtractDCoffset =0;  % subtract DC OFFSET
        params.removeEYE        =1;  % remove EYE artifacts
    case 7
        params.NOTCHfilter      =0;  % apply NOTCH filter
        params.LPfilter         =0;  % apply   LP  filter
        params.SubtractMean     =0;  % subtract mean values
        params.SubtractDCoffset =1;  % subtract DC OFFSET
        params.removeEYE        =0;  % remove EYE artifacts
end


end
