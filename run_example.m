clear;
close all;
 

load('iML1515.mat')
model=iML1515;

% get biomass reaction
model.csense=columnVector(model.csense);
biomassRxn=model.rxns{model.c==1};


% define target produc
%targetRxn='EX_succ_e';
targetRxn='EX_lyco_e';
targetID=findRxnIDs(model,targetRxn);

% set uptake rate of oxygen and glucose
oxygenRxn='EX_o2_e';
substrate='EX_glc__D_e';
% model = changeRxnBounds(model,{substrate,oxygenRxn},[-20,-20],{'l','l'});

% limit reaction rate in realistic range
model.lb(model.lb<-100)=-100;
model.ub(model.ub>100)=100;

orimodel=model;

% when the model size is big, it is better to compress the model so that
% the linear reactions can be reduced
if size(model.S,1)>200
    % compress the model and get compressed candidate reactions for knockout
    [model,candidate]=nihba_prep(orimodel,substrate,oxygenRxn,biomassRxn,targetRxn);
end
%candidate.rxns=model.rxns;

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
options.bendersTime=3600*60; % two hours

% minimum ratio of growth in production mutants
minRationOfGrowth=0.1; % 10% wild type growth

% call run_nihba and save results in excel and mat files
run_nihba(model,options, minRationOfGrowth);