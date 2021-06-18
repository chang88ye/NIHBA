clear;
close all;
% 
% load('e_coli_core.mat');
% model=e_coli_core;

% load('iAF1260b.mat');
% model=iAF1260b;

load('iML1515.mat')
%model=iML1515;

% load('narg_ecoli_model_pw0.mat')
% model=modelOut;

model.csense=model.csense';
biomassRxn=model.rxns{model.c==1};

targetRxn='EX_succ_e';
targetID=findRxnIDs(model,targetRxn);

% tilting strategy used in inner LP

oxygenRxn='EX_o2_e';
substrate='EX_glc__D_e';
% substrate='EX_glyc_e';
% substrate='EX_xyl__D_e';
% model = changeRxnBounds(model,'EX_glc__D_e',0,'l');
% model = changeRxnBounds(model,{substrate,oxygenRxn},[-20,-20],{'l','l'});

% model.lb(model.lb<-10000)=-10000;
% model.ub(model.ub>10000)=10000;
orimodel=model;

[model,candidate]=nihba_prep(orimodel,substrate,oxygenRxn,biomassRxn,targetRxn);
 
model.lb(model.lb<-100)=-100;
model.ub(model.ub>100)=100;
optimizeCbModel(model)

candidateTargets={'EX_mal__L_e/MALtex', 'EX_succ_e','EX_etoh_e/ETOHtex/ALCD2x/ETOHtrpp','EX_ac_e/ACtex','EX_fum_e/FUMtex'};%'EX_mal__L_e/MALtex','EX_ala__L_e/ALAtex'};

targetRxn=candidateTargets{2};
selectedRxns=setdiff(candidate.rxns, {'ATPM', biomassRxn, targetRxn});
disp(['The size of candidates for knockout is: ', num2str(length(selectedRxns))]);


geneDelFlag = false;
nPts = 30;
figure(1)
colors = jet(25);
hold on
solWT=optimizeCbModel(model);
maxGrowth=solWT.f;
productionEnvelope(model,{},'r',targetRxn,biomassRxn,geneDelFlag,nPts);
xlabel('Biomass', 'FontSize', 20);
ylabel('Production Rate', 'FontSize', 20);

values=[];
koSol={};

options.selectedRxns=selectedRxns;
options.targetRxn=targetRxn;
options.biomassRxn=biomassRxn;
options.innerRxns=targetRxn;
options.innerOsense='max';
options.maxKO=5;
options.minKO=2;
options.timeLimit=3600*2;
options.bendersTime=3600*2;

tmodel = changeRxnBounds(model,biomassRxn,0.1*maxGrowth,'l');

[solutions,~]=run_nihba(tmodel, options);


productionEnvelope(model,solutions.koSet, 'b',targetRxn,biomassRxn,geneDelFlag,nPts);


for i=1:length(solutions.allSet)
    %for i=1:21
    deletions=solutions.allSet{i};
    [~, maxGrowth, maxProd, minProd]=analyzeOptKnock(model,deletions, targetRxn);
    values(end+1,:)=[maxGrowth, minProd, maxProd];
    koSol(end+1,:)=[deletions',cell(1,options.maxKO-length(deletions))];
    hold on
    productionEnvelope(model,deletions,'b',targetRxn,biomassRxn,geneDelFlag,nPts);
end

tabnums=array2table(values,'VariableNames',{'biomass','minProd','maxProd'});
tabstrs=cell2table(koSol, 'VariableNames',cellfun(@(x) ['ko' num2str(x)],num2cell(1:options.maxKO),'UniformOutput',false));
tabs=[tabnums,tabstrs];

tmpSet=strsplit(targetRxn,'/');
filename=[tmpSet{1}, '_KO', num2str(options.maxKO), '_', model.description,...
    '_',solutions.method,'-',num2str(length(selectedRxns))];
saveResults(tabs, filename);

if isfield(solutions, 'EC')
    EC=solutions.EC;
    save([filename,'_EC.mat'],'EC');
end


%% subfunctions

function saveResults(tabVal, filename)
if ~exist('results', 'dir')
    mkdir('results');
end

writetable(tabVal,['results', filesep, filename,'.csv']);
save(['results', filesep, filename, '.mat'], 'tabVal');
end
