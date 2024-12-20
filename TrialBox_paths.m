function   paths=TrialBox_paths(code_dir)
% function paths=TrialBox_paths(data_dir)
    
p='';
d=pathsep;%';';
S=filesep;

if nargin < 1
    code_dir = pwd;
end
p=[p code_dir 'main'               d]; %% main (executables)
% p=[p code_dir 'DL'                 d]; %% dictionary learning tests (executables)
p=[p genpath([ code_dir 'DL'     ] )]; %% dictionary learning tests (executables)
p=[p genpath([ code_dir 'models' ] )]; %% models load
% p=[p code_dir 'models'             d]; 
% p=[p code_dir 'models'    S 'pms'  d]; %% model params load
p=[p genpath([ code_dir 'filters' ] )]; %% models load
% p=[p code_dir 'filters'            d]; %% filters load
% p=[p code_dir 'filters'   S 'pms'  d]; %% filter params load
p=[p code_dir 'plots'              d]; %% plots load
p=[p code_dir 'plots'     S 'pms'  d]; %% plot params load
p=[p code_dir 'utilities'          d]; %% plots load
p=[p code_dir 'utilities' S 'pms'  d]; %% plots load
p=[p code_dir 'datasets'           d]; %% data load
p=[p code_dir 'datasets'  S 'Graz' d]; %% Graz data load
addpath(p);
if nargout>0
    paths=p;
end
