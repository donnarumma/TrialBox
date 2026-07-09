function [Pos, Val]=Odysseas_Get_Nearest_Value(Number,InArray)
%Gets the nearest value of an array to a number
Pos=NaN(size(Number,1),1);
Val=NaN(size(Number,1),1);

for i=1:size(Number,1)
    [~,Pos(i)]=min(abs(InArray-Number(i)));
    Val(i)=InArray(Pos(i));
end
end