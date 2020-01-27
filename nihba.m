function solutions=nihba(model, options)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% NIHBA: A NETWORK INTERDICTION APPROACH WITH HYBRID BENDERS ALGORITHM FOR STRAIN DESIGN
% solve:
%       min D*y
%       s.t. A*x+B*y >=b; x in {0, 1}, y>=0;

% USAGE:
%
%    solutions=NIHBA(model, options)
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
%                      
% OUTPUTS:
%    solution:               normal solution structure, which additionally includes:
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
MAXFLUX = 100;
MAXDUAL = 100;

if (~exist('options','var') || isempty(options) )
    error('No target reaction specified')
else
    %Default Values
    if isfield(options,'targetRxn')
        selTargetRxns = logical(ismember(model.rxns,options.targetRxn));
        if ~any(selTargetRxns)
            error([options.targetRxns ' not found. Double check spelling.']);
        end
    else
        error('No target reaction specified');
    end
    
    if ~isfield(options,'x0'), options.x0=[]; end
    if ~isfield(options,'preSol'), options.preSol=[]; end
    if ~isfield(options,'method'), options.method='sum'; end
    if ~isfield(options,'innerRxns'), options.innerRxns=model.rxns(model.c==1); end
    if ~isfield(options,'innerOsense'), options.innerOsense='max'; end
    if ~isfield(options,'selectedRxns'), options.selectedRxns=model.rxns;end
    if ~isfield(options,'maxKO'), options.maxKO=5; end
    if ~isfield(options,'minKO'), options.minKO=2; end
    if ~isfield(options,'koType'), options.koType='rxns'; end
    
    % import options to trigger truncated branch and bound
    if ~isfield(options,'iterationLimit'), options.iterationLimit=1e30; end
    if ~isfield(options,'timeLimit'), options.timeLimit=300 ; end
    if ~isfield(options,'Heuristics'), options.Heuristics = 1.0; end
    if ~isfield(options,'MIPFocus'), options.MIPFocus = 1; end
    if ~isfield(options,'ImproveStartGap'), options.ImproveStartGap = Inf; end
    
    
    % options for Benders algorithm
    if ~isfield(options,'bendersTime'), options.bendersTime=1e30;end
    if ~isfield(options,'bendersIteration'), options.bendersIteration=1e30;end
end

% solver options
solverOpt.iterationLimit=options.iterationLimit;
solverOpt.timeLimit=options.timeLimit;
solverOpt.Heuristics=options.Heuristics;
solverOpt.MIPFocus=options.MIPFocus;
solverOpt.ImproveStartGap=options.ImproveStartGap;

if isfield(options,'innerRxns')
    selInnerRxns = logical(ismember(model.rxns,options.innerRxns));
    if ~any(selInnerRxns)
        error([options.innerRxns ' not found. Double check spelling.']);
    end
else
    warning('No target specified for inner optimisation, minisation of sum of fluxes is assumed.')
    selInnerRxns = model.c;
end

if isfield(options,'selectedRxns')
    selSelectedRxns = logical(ismember(model.rxns,options.selectedRxns));
end

switch lower(options.koType)
    case 'rxns'
        %% Generate selection reaction matrix
        model.selRxnMatrix = selMatrix(selSelectedRxns)';
        possibleKOList = model.rxns(selSelectedRxns);
        
    case 'genesets'
        %% Generate reaction gene set mapping matrix
        %remove biomass reaction from grRules and generate unique gene set list
        possibleKOList = unique(model.grRules(selSelectedRxns));
        if isempty(possibleKOList{1}), possibleKOList = possibleKOList(2:end); end
        for i = 1:length(possibleKOList)
            model.selRxnMatrix(:,i) =  double(strcmp(possibleKOList{i},model.grRules));
        end
        
    case 'genes'
        %% Use rxnGeneMat as selRxnMatrix
        model.selRxnMatrix = model.rxnGeneMat;
        possibleKOList = model.genes;
    otherwise
        error('Unrecognized KO type')
end

%% Generate koCost if not present
if ~isfield(options,'koCost')
    if isfield(model,'koCost')
        if length(model.koCost) == 1
            options.koCost = ones(length(possibleKOList,1)) * model.koCost;
        else
            options.koCost = model.koCost;
        end
    else
        options.koCost = ones(length(possibleKOList),1);
    end
elseif length(model.koCost) == 1
    options.koCost = ones(length(possibleKOList,1)) * model.koCost;
else
    options.koCost = model.koCost;
end

%% Setup model
model.ub(isinf(model.ub)) = MAXFLUX;
model.ub(model.ub>MAXFLUX) = MAXFLUX;
model.lb(isinf(model.lb)) = -MAXFLUX;
model.lb(model.lb<-MAXFLUX) = -MAXFLUX;

% add pertubation to unbounded flux [0.9,1.1]
maxFlux=max(model.ub);
tmpBool=(model.ub==maxFlux);
model.ub(tmpBool)=maxFlux*(0.9+0.2*rand([sum(tmpBool),1]));
tmpBool=(model.lb==-maxFlux);
model.lb(tmpBool)=-maxFlux*(0.9+0.2*rand([sum(tmpBool),1]));

maxDual=MAXDUAL*(0.9+0.2*rand([length(model.rxns),1]));

% % maxDual=MAXDUAL;


model.rxnsPresent = ones(length(model.rxns),1);

tmodel=changeObjective(model,options.targetRxn);
tsol=optimizeCbModel(tmodel);
TMP=tsol.f;
%% Create bi-level MILP problem
[nMets, nRxns] = size(model.S);
nInt = size(model.selRxnMatrix,2);

% if strcmp(options.innerRxns, options.biomassRxn)
%     selInnerRxns=-selInnerRxns;
% % else
% %     selInnerRxns=1e-3*selInnerRxns-model.c;
% end

% Initialisation
epsilon=1e-4;


[x0,xc]=generateCorePoint(model, selTargetRxns, options);
% 
% x0=zeros(length(xc),1);
% x0([94, 159, 168, 249, 437])=1;
scentre=x0;


lambda=0.5;% parameters
mu=1e-8;

k=2;
nc=2*options.maxKO;

xBest=[];
xSet=[];

LB=-1e30;
UB=1e30;
zu=UB;
zl=LB;

subObj=@(x)[zeros(nMets,1);
    double(selInnerRxns);
    -model.lb.*(1-model.selRxnMatrix * x);
    model.ub.*(1-model.selRxnMatrix * x);
    maxDual.*(model.selRxnMatrix*x);
    maxDual.*(model.selRxnMatrix*x);
    0;
    ];

% solutions creating feasibility cuts
Fset=[];
% dual vectors of Fset
Fdual=[];

% solutions creating optimality cuts
% Oset=[];
% dual vectors of Oset
Ocuts=struct();
Ocuts.A=[];
Ocuts.b=[];
Ocuts.csense=char();

Fcuts=struct();
Fcuts.A=[];
Fcuts.b=[];
Fcuts.csense=char();

% inform gurobi to tell infeasible or unbounded
params.InfUnbdInfo = 1;
params.Method=0; % prime simplex is the fastest in our problem
subBasic=createSubproblem(model, selTargetRxns, selInnerRxns);
subBasic.c=zeros(nMets+3*nRxns+2*nRxns+1,1);
sol=solveCobraLP(subBasic);
subBasic.basis=sol.basis;

SMP0=createStablisedMasterProblem(model,options);

ii=0;
iteration=1;
optimality=false;
solveSP=true;
mipStart=[]; % start solution for mip

iEC=[]; % iterational evolutionary curve
tEC=[]; % time evolutionary curve

tic;
while (~optimality) && (toc<=options.bendersTime)
    
    if solveSP
        % update subproblem
        c1=subObj(x0);
        c2=subObj(xc);
        
        % update objective of dual subproblem
        subBasic.c=c1+mu*c2;
        
        % solve the dual subproblem
        sol=solveCobraLP(subBasic, params);
  
        % update master constraints
        if sol.stat==1 % bounded
            xc=lambda*(xc+x0);
            
            zl=c1'*sol.full;
            
            if LB<zl
                LB=zl;
                xBest=x0;
            end
            
            if zl>0.1 % threshold
                xSet=[xSet; x0'];
            end
            
            if  zu-zl<epsilon || zu<LB || zu<0.1*TMP
                sol=solveCobraMILP(SMP);
                if abs(sol.obj-LB)<epsilon % output optimal LB
                    optimality=true;
                end
                % reverse local branching constraint
                SMP0.A=[
                    SMP0.A;
                    ((scentre == 0) - (scentre == 1))' 0; % logical cut
                    ];
                SMP0.b(end+1)=k+1-sum(scentre==1);
                SMP0.csense(end+1)=char('>');
                
                scentre=x0;
                % reset k;
                %k=1;
            else
                u=sol.full(nMets+1:nMets+nRxns);
                s=sol.full(nMets+nRxns+1:nMets+2*nRxns);
                t=sol.full(nMets+2*nRxns+1:nMets+3*nRxns);
                r1=sol.full(nMets+3*nRxns+1:nMets+4*nRxns);
                r2=sol.full(nMets+4*nRxns+1:nMets+5*nRxns);
                Ocuts.A=[
                    Ocuts.A;
                    (model.ub.*t -model.lb.*s-maxDual.*r1-maxDual.*r2)'*model.selRxnMatrix 1; % optimality cut
                    ];
                
                Ocuts.b(end+1,:)=-model.lb'*s+model.ub'*t + selInnerRxns'*u;
                Ocuts.csense(end+1,:)='<';
            end
            
        else %sol.stat==2 % unbounded
            try
                u=sol.full(nMets+1:nMets+nRxns);
                s=sol.full(nMets+nRxns+1:nMets+2*nRxns);
                t=sol.full(nMets+2*nRxns+1:nMets+3*nRxns);
                r1=sol.full(nMets+3*nRxns+1:nMets+4*nRxns);
                r2=sol.full(nMets+4*nRxns+1:nMets+5*nRxns);
                
                Fset(end+1,:)=x0;
                Fdual(end+1,:)=sol.full;
                
                ii=ii+1;
                Fcuts.A=[
                    Fcuts.A;
                    (model.ub.*t -model.lb.*s-maxDual.*r1-maxDual.*r2)'*model.selRxnMatrix 0; % feasibility cut
                    ];
                
                Fcuts.b(end+1,:)=-model.lb'*s+model.ub'*t + selInnerRxns'*u;
                Fcuts.csense(end+1,:)='<';
            catch
                params.NumericFocus=3;
                params.Method=-1;
                disp('Warning: no extreme ray returned. Please set the gurobi output as x=resultgurobi.unbdray.');
            end
        end
    end
%     if mod(iteration,100)==0
%         [Fset,Fdual]=findMinimalFeasibilityCut(subObj, Fset, Fdual, Fcuts);
%     end

    
    SMP=updateMasterProblem(SMP0, Ocuts,Fcuts);
    
    SMP1=SMP;
    %%local branching constraint
    SMP1.A=[
        SMP1.A;
        ((scentre == 0) - (scentre == 1))' 0; % logical cut
        ];
    SMP1.b(end+1)=k-sum(scentre==1);
    SMP1.csense(end+1)=char('<');
    SMP1.x0=mipStart;
    
    solverOpt.MIPGap=300/(iteration.^0.5+1)+1;
    sol=solveCobraMILP(SMP1,solverOpt);

    if (sol.stat==0)%infeasible
        if k>=nc
            optimality=true;
            break;
        end
        
        % reverse local branching constraint
        SMP0.A=[
            SMP0.A;
            ((scentre == 0) - (scentre == 1))' 0; % logical cut
            ];
        SMP0.b(end+1)=k+1-sum(scentre==1);
        SMP0.csense(end+1)=char('>');
        % reset k;
        k=k+1;
        solveSP=false;
    else %(sol.stat==1 || ~isempty(sol.obj)) && (~all((sol.int>1e-4)==scentre)) && (sol.cont>=LB)%(sol.cont>0)
        x0=sol.int>1e-4;
        zu=sol.cont;
        solveSP=true;
        mipStart=sol.full;
    end
    iEC(end+1)=LB;
    tEC(end+1)=toc;
    
    if iteration==1 || mod(iteration,100)==0
        disp([num2str(iteration), ' iteration(s) passed, cpu time = ', num2str(toc),' seconds.']);
        disp(['LB= ', num2str(LB), ' ,UB= ', num2str(zu)]);
    end
    iteration=iteration+1;
end
solutions.x=xBest;
solutions.obj=LB;
solutions.koSet=possibleKOList(logical(xBest));
nsol=size(xSet,1);
for i=1:nsol
    solutions.allSet{i}=possibleKOList(logical(xSet(i,:)));
end
solutions.allSet{nsol+1}=solutions.koSet;

solutions.allSet=solutions.allSet(~cellfun('isempty',solutions.allSet));
solutions.EC{1}=iEC;
solutions.EC{2}=tEC;
solutions.method='NIHBA';
end


function LPproblem=createStablisedMasterProblem(model,options)
% [nMets, nRxns] = size(model.S);
nInt=size(model.selRxnMatrix,2);
% [y, z]
LPproblem.A=[
    ones(1,nInt) 0;
    ones(1,nInt) 0;
    %     H*ones(1,nInt) -1;
    ];
LPproblem.b=[
    options.maxKO;
    options.minKO
    %     H;
    ];
LPproblem.csense=char([
    '<';
    '>';
    %     '>'
    ]);

LPproblem.c=[zeros(nInt, 1); 1];

LPproblem.lb(1:nInt)=0;
LPproblem.ub(1:nInt)=1;
LPproblem.lb(end+1)=0;
LPproblem.ub(end+1)=Inf;

LPproblem.vartype=char([
    'B'*ones(nInt,1);
    'C';
    ]);
LPproblem.osense=-1; % maximisation
LPproblem.x0=[];
end


function subBasic=createSubproblem(model, selTargetRxns, selInnerRxns)
[nMets, nRxns] = size(model.S);
H=Inf;

innerSign=1;
if any(selInnerRxns<0) % indicating maximisation for inner problem
    innerSign=-1;
end

% create subproblem
% [k u, s, t, r1, r2, beta]
subBasic.A=[
    sparse(nMets,nMets)     model.S                 sparse(nMets, 2*nRxns+1+2*nRxns);
    model.S'                sparse(nRxns,nRxns)     -speye(nRxns,nRxns)          speye(nRxns,nRxns+2*nRxns)     double(selInnerRxns);
    sparse(nRxns,nMets)     speye(nRxns,nRxns)     sparse(nRxns, 2*nRxns+2*nRxns)         -model.lb;
    sparse(nRxns,nMets)     -speye(nRxns,nRxns)      sparse(nRxns, 2*nRxns+2*nRxns)        model.ub;
    sparse(nRxns,nMets)     innerSign*speye(nRxns,nRxns)      sparse(nRxns, 2*nRxns)       speye(nRxns,nRxns)  -speye(nRxns,nRxns+1);
    ];
subBasic.b=[
    model.b;
    double(selTargetRxns);
    zeros(nRxns, 1);
    zeros(nRxns, 1);
    zeros(nRxns, 1);
    ];
subBasic.csense=char([
    '=' * ones(nMets, 1);
    '=' * ones(nRxns, 1);
    '>' * ones(nRxns, 1);
    '>' * ones(nRxns, 1);
    '=' * ones(nRxns, 1);
    ]);

subBasic.lb=[
    -H * ones(nMets, 1);
    -H * ones(nRxns, 1);
    zeros(nRxns,1);
    zeros(nRxns,1);
    zeros(nRxns,1);
    zeros(nRxns,1);
    0;
    ];
subBasic.ub=[
    H * ones(nMets, 1);
    H * ones(nRxns, 1);
    H * ones(nRxns, 1);
    H * ones(nRxns, 1);
    H * ones(nRxns,1);
    H * ones(nRxns,1);
    H;
    ];
subBasic.vartype=char('C'*ones(nMets+3*nRxns+1+2*nRxns,1));
subBasic.osense=1; % minimisation
subBasic.x0=[];
end


function [x0, corepoint]=generateCorePoint(model, selTargetRxns, options)
[nMets, nRxns] = size(model.S);
nInt=size(model.selRxnMatrix,2);

% create subproblem
% [v, y]
LPproblem.A=[
    model.S                 sparse(nMets,nInt);
    speye(nRxns)            model.selRxnMatrix.* repmat(model.ub, 1, nInt);
    speye(nRxns)            model.selRxnMatrix.* repmat(model.lb, 1, nInt);
    sparse(1,nRxns)         ones(1,nInt);
    sparse(1,nRxns)         ones(1,nInt);
    ];
LPproblem.b=[
    model.b;
    model.ub;
    model.lb;
    options.maxKO;
    options.minKO
    ];
LPproblem.csense=char([
    '='*ones(nMets,1);
    '<'*ones(nRxns,1);
    '>'*ones(nRxns,1);
    '<';
    '>';
    ]);

LPproblem.lb(1:nRxns)=model.lb;
LPproblem.ub(1:nRxns)=model.ub;
LPproblem.lb(nRxns+1:nRxns+nInt)=0;
LPproblem.ub(nRxns+1:nRxns+nInt)=1;

LPproblem.vartype=char([
    'C'*ones(nRxns,1);
    'B'*ones(nInt,1);
    ]);
LPproblem.osense=-1; % maximisation
LPproblem.x0=[];

LPproblem.c=[selTargetRxns; zeros(nInt, 1); ];

solutionPool=[];
solution=solveCobraMILP(LPproblem);
y0=solution.int;
solutionPool(:,1)=y0;
x0=y0;

cLPproblem=LPproblem;
for c=1:5 % get some feasible solutions
    cLPproblem.A(end+1,:)=[sparse(1,nRxns) ((y0 == 0) - (y0 == 1))' ]; % logical cut]
    cLPproblem.b(end+1)=1-sum(y0==1);
    cLPproblem.csense(end+1)=char('>');
    
    solution=solveCobraMILP(LPproblem);
    y0=solution.int;
    solutionPool(:,end+1)=y0;
end
corepoint=mean(solutionPool,2);
end


function masterProblem=updateMasterProblem(masterProblem, Ocuts,Fcuts)

masterProblem.A=[
    masterProblem.A;
    Ocuts.A;
    Fcuts.A;
    ];


masterProblem.b=[
    masterProblem.b;
    Ocuts.b;
    Fcuts.b;
    ];

masterProblem.csense=[
    masterProblem.csense;
    Ocuts.csense;
    Fcuts.csense;
    ];
end