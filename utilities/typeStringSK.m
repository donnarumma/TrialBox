function  str=typeStringSK(tType1,tType2)
%function str=typeStringSK(tType1,tType2)

if     tType1==1
  s='ACT S, OBS K';     
elseif tType1==2
  s='OBS S, ACT K';      
elseif tType1==3
  s='ACT S, ACT K';       % Joint Action
% elseif tType1==4
%   s='EYE';      % Eye
% elseif tType1==5
%   s='CONTROL';  % Control
end  

str=sprintf('%s (D%g)',s,tType2);