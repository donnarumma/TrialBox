function   X = Graz2A_getData(data,time,tmin,tmax,fc)
% extraction data for class
% function X = Graz2A_getData(data,time,tmin,tmax,fc)
    X           = cell(length(time),1);
    data        = data';
    it_start    = time+tmin*fc;
    it_stop     = time+tmax*fc-1;
    if it_stop < size(data,2)
        X{1,1}  = data(1:22,it_start:it_stop);
    else
        X{1,1}  = data(1:22,it_start:size(data,2));
    end
end