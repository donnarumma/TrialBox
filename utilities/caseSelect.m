function [it1,it2] = caseSelect(c)

t_s4s = [0.8998	1.5429	1.3018	1.0606	1.4626	1.2214	1.3822	1.3018	1.5429];

t_2cl = [0.9952	1.5468	1.3047	1.1488	1.5508	0.9717	1.3004	1.8104	1.5468];

t_scl = [0.8998	1.5429	1.3018	1.1410	1.5429	0.9802	1.3018	1.8645	0.7391];


if c==1
    % case 1 
    % tsub1 = T*s4s
    % tsub3 = T*2cl
    it1 = t_s4s;
    it2 = t_2cl;
elseif c==2
    % case 2 
    % tsub1 = 0.5-2.5
    % tsub3 = T*scl
    it1 = 0;
    it2 = t_scl;
elseif c==3
    % case 3
    % tsub1 = 0.5-2.5
    % tsub3 = T*s4s
    it1 = 0;
    it2 = t_s4s;
elseif c==4
    % case 4 
    % tsub1 = 0.5-2.5
    % tsub3 = T*2cl
    it1 = 0;
    it2 = t_2cl;
end