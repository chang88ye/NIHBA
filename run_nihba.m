function [solutions, bestKO]=run_nihba(model,options, minRatioOfGrowth)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% NIHBA: A NETWORK INTERDICTION APPROACH WITH HYBRID BENDERS ALGORITHM FOR STRAIN DESIGN
% solve:
%       min D*y
%       s.t. A*x+B*y >=b; x in {0, 1}, y>=0;

% USAGE:
%
%    run_nihba(model,options, minRationOfGrowth)
%
%
% INPUTS:
%    model:                  Cobra model structure
%    options:                a parameter structure, members include:
%                              *  substrateRxn - carbon source uptake reaction ('EX_glc__D_e' by default)
%                              *  selectedRxns - a list of candidate reactions that can be knocked out
%                              *  targetRxn    - target reaction for which flux increase leads to high chemical production
%                              *  biomassRxn   - biomass reaction 
%                              *  innerRxns    - reaction to be optimised in the inner level of the bilevel framework 
%                              *  minKO        - minimum number of knockable reactions (default: 0)
%                              *  maxKO        - maximum number of knockable reactions (default: 5) 
%                              *  timeLimit    - run time in seconds (default: 252000)
%                              *  x0           - initialisation of knockout reactions
%                              *  preSol       - solution set that has been found but required not to appear in new solutions
%   minRationOfGrowth:       minimum ratio of growth required for survival
%                      
% OUTPUTS: save results in mat file, a  normal solution structure, which additionally includes:
%                               * koSet         - best knockout set returned 
%                               * allSet        - all sets of knockout
%                               * EC            - evolutionary trajectory
%                               * method        - method used
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Authors:
%       - Shouyong Jiang 31/01/19
%       - inquires about the implimentation should be directed to math4neu@gmail.com
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if nargin<3
    % set minimum ratio of growth to 0.1 for production mutants
    minRatioOfGrowth=0.1;  
end

% set requirements for production mutants
solWT=optimizeCbModel(model);
maxGrowth=solWT.f;
tmodel = changeRxnBounds(model,model.rxns{model.c==1},minRatioOfGrowth*maxGrowth,'l');

% apply nihba to find manipulation strategies
solutions=nihba(tmodel, options);

% calculate each mutant's phenotype
values=[];
koSol={};
for i=1:length(solutions.allSet)
    deletions=solutions.allSet{i};
    [~, maxGrowth, maxProd, minProd]=analyzeOptKnock(model,deletions, options.targetRxn);
    values(end+1,:)=[maxGrowth, minProd, maxProd];
    koSol(end+1,:)=[deletions',cell(1,options.maxKO-length(deletions))];
end

% find best strategy in terms of lower bound production
[~,idx]=max(values(:,2));
tmpBestKO=solutions.allSet{idx};
bestKO={};
for i=1:length(tmpBestKO)
    tmpSet=strsplit(tmpBestKO{i},'/');
    bestKO{i}=tmpSet{1};
end

% creat a table of results
tabnums=array2table(values,'VariableNames',{'biomass','minProd','maxProd'});
tabstrs=cell2table(koSol, 'VariableNames',cellfun(@(x) ['ko' num2str(x)],num2cell(1:options.maxKO),'UniformOutput',false));
tabs=[tabnums,tabstrs];

% save the table
tmpSet=strsplit(options.targetRxn,'/');
filename=[tmpSet{1}, '_KO', num2str(options.maxKO), '_', model.description,...
    '_',solutions.method,'-',num2str(length(options.selectedRxns))];
saveResults(tabs, filename);

% save the runtime
if isfield(solutions, 'EC')
    EC=solutions.EC;
    save([filename,'_EC.mat'],'EC');
end
end

%% subfunctions
function saveResults(tabVal, filename)
if ~exist([pwd,'results'], 'dir')
    mkdir('results');
end

writetable(tabVal,['results', filesep, filename,'.csv']);
save(['results', filesep, filename, '.mat'], 'tabVal');
end
