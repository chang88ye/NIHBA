function nihba_setup()

% launch cobratoolbox
try 
    initCobraToolbox
catch
    error(' > cobratoolbox is not available.\n');
end

% use Gurobi (if installed) as the default solver for LP, QP and MILP problems
solOK=changeCobraSolver('gurobi', 'ALL', 1);
if ~solOK
    fprintf(2,'***GUROBI is not fully accessible, please check gurobi installation or switch to another MILP solver.\n')
end

% make sure the current folder is NIHBA
str = which('nihba.m');
if isempty(str)
    fprintf(2,'***Current directory is not NIHBA, please change it to NIHBA.\n')
end

% add NIHBA to matlab search path
curPath=fileparts(str);
addpath(curPath);
end