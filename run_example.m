clear;
close all;

% use E. coli core model
load('ecoli_core_model.mat');

% use latest E. Coli Model
%load('iML1515.mat')

% get biomass reaction
%model.csense=columnVector(model.csense);
biomassRxn=model.rxns{model.c==1};


% define target produc
targetRxn='EX_succ_e';
%targetRxn='EX_lyco_e';

% set uptake rate of oxygen and glucose
oxygenRxn='EX_o2_e';
substrate='EX_glc__D_e';
% model = changeRxnBounds(model,{substrate,oxygenRxn},[-20,-20],{'l','l'});

if strcmp(lower(model.description),'ecoli_core_model') % identifier is different in ecoli_core model
    targetRxn='EX_succ(e)';
    oxygenRxn='EX_o2(e)';
    substrate='EX_glc(e)';
end


% limit reaction rate in realistic range
model.lb(model.lb<-100)=-100;
model.ub(model.ub>100)=100;

orimodel=model;

% when the model size is big, it is better to compress the model so that
% the linear reactions can be reduced
if size(model.S,1)>200
    % compress the model and get compressed candidate reactions for knockout
    [model,candidate]=nihba_prep(orimodel,substrate,oxygenRxn,biomassRxn,targetRxn);
else
    candidate.rxns=model.rxns; 
end

% make sure target reaction is in the reaction list of compressed model
if ~strcmp(model.rxns, targetRxn)
    targetRxn=model.rxns{contains(model.rxns, targetRxn)};
end

% remove unreasonable reactions from candidate knockout set
selectedRxns=setdiff(candidate.rxns, {'ATPM', biomassRxn, targetRxn});

disp(['Preprocessing step completed, ', num2str(length(selectedRxns)),' candidate reactions are selected for knockout.'])

% sett options for nihba algorithm
options.selectedRxns=selectedRxns;
options.targetRxn=targetRxn;
options.innerRxns=targetRxn;
options.innerOsense='max'; %-------------
options.maxKO=5; % at most the number of knockouts
options.minKO=2; % at least two knockouts 
options.timeLimit=3600*2;   % two hours
options.bendersTime=30; %3600*60; % two hours

% minimum ratio of growth in production mutants
minRationOfGrowth=0.1; % 10% wild type growth

% call run_nihba and save results in excel and mat files
[~,best_strategy]=run_nihba(model,options, minRationOfGrowth);

% Draw the flux values on the map "target.svg" which can be opened in FireFox
if strcmp(lower(orimodel.description),'ecoli_core_model')
    map=readCbMap('ecoli_Textbook_ExportMap');
    %options.lb = -10;
    %options.ub = 10;
    options.zeroFluxWidth = 0.1;
    options.rxnDirMultiplier = 10;
    options.fileName='WildType.svg';
    %options.textsize = 12;
    
    % wild type
    WTsolution=optimizeCbModel(orimodel);
    drawFlux(map, orimodel, WTsolution.x, options);
    web(options.fileName);
    
    % mutant strain
    tmodel=changeRxnBounds(orimodel,best_strategy,0);
    MTsolution=optimizeCbModel(tmodel);
    options.fileName='Mutant.svg';
    drawFlux(map, orimodel, MTsolution.x, options);
    web(options.fileName, '-new');
end