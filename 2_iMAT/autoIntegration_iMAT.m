function [PFD,N_highFit,N_zeroFit,minLow,minTotal,OpenedLReaction,wasteDW,Hreactions,Lreactions,MILP] = autoIntegration_iMAT(model,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCatag,doMinPFD,type)
%% inputs
% model ... the model; model needs to be constrained
% storeProp,SideProp ... proportion of sides and storage molecules
% epsilon_f/r ... calculated epsilon vector for every reaction
% ExpCatag ... the expressiom category of each gene; note that the gene
% name should in iCEL format

%apply default
if (nargin < 8)
    doMinPFD = true;
end
if (nargin < 9) %fullsensitivity means sensitivity that includes negative sensitivity score (confident no flux), and due to technical reason, in full sensitivity, -1/+1 score means violation of minTotal flux
    type = 'canonical';
end
%set global constant 
bacMW=966.28583751;
%changeCobraSolverParams('LP','optTol', 10e-9);
%changeCobraSolverParams('LP','feasTol', 10e-9);
tol = 1e-5; %it could be 1%
%%
fprintf('Start flux fitting... \n');
tic()
fprintf('Processing expression data... \n');
worm = model;
% process gene expression data
expression_gene=struct;
expression_gene.value = [3*ones(length(ExpCatag.high),1);2*ones(length(ExpCatag.dynamic),1);1*ones(length(ExpCatag.low),1);zeros(length(ExpCatag.zero),1)];
expression_gene.gene = [ExpCatag.high;ExpCatag.dynamic;ExpCatag.low;ExpCatag.zero];
[expressionRxns] = mapExpressionToReactions_xl(worm, expression_gene);
%% Step1: do iMAT integration
toc()
fprintf('fitting iMAT... \n');
tic()
% adjust the ATPm according to the bacterial uptake. This requires the
% manual tuning
worm = changeRxnBounds(worm,'RCC0005_I',ATPm,'l');
worm = changeRxnBounds(worm,'RCC0005_X',ATPm,'l');
% change the side proportion
worm.S(end-1, strcmp('EXC0050_L',worm.rxns)) = storeProp*bacMW*0.01;
worm.S(end, strcmp('EXC0050_L',worm.rxns)) = SideProp*bacMW*0.01;
if strcmp(type,'canonical') %low category is merged with zero (rare) so that there is NO minimization of low flux
    [solution, MILProblem,Hreactions,Lreactions] = iMAT_xl(worm, expressionRxns, epsilon_f, epsilon_r,1.1,2.9);
elseif strcmp(type,'iMATplus')
    [solution, MILProblem,Hreactions,Lreactions] = iMAT_xl(worm, expressionRxns, epsilon_f, epsilon_r,0.1,2.9);
else
    error('invalid integration type!');
end
toc()
% Step2: do minimization of low reactions
%set solution as constraint
if strcmp(type,'iMATplus')
    fprintf('Doing the iMAT_plus, minimizing low flux... \n');
    tic()
    c_ori = zeros(size(MILProblem.A,2),1); %log down the original variable number
    MILProblem = solution2constraint(MILProblem,solution);
    MILProblem = addAbsFluxVariables(MILProblem, worm);
    %minimize low flux 
    %to save computation labor, we do the flux minimization by V+ we generated
    %in step 1 matrix 
    %create a new objective function
    % Creating c (objective function)
    c_minFluxLowRxns = zeros(length(worm.rxns),1);
    c_minFluxLowRxns(expressionRxns == 1 | expressionRxns == 0) = 1; %set the coeffi of low and zero reactions as 1
    c = [c_ori;c_minFluxLowRxns];
    MILProblem.c = c;
    % set the objective sense as minimize 
    MILProblem.osense = 1;
    % sovle it
    solution = solveCobraMILP_XL(MILProblem, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 0);
    if solution.stat ~= 1
        error('infeasible or violation occured!');
    end
    toc()
    fprintf('Minimizing low flux completed! \n');
    minLow = solution.obj;
    solution.obj = solution.obj + tol;
    MILProblem = solution2constraint(MILProblem,solution);
else
    MILProblem = solution2constraint(MILProblem,solution);
    c_ori = zeros(size(MILProblem.A,2),1); %log down the original variable number
    MILProblem = solution2constraint(MILProblem,solution);
    MILProblem = addAbsFluxVariables(MILProblem, worm);
    minLow = 0;
end
% step3: minimize total flux as needed
MILP = MILProblem; %write the MILP output, this MILP contains all data-based constraints (minTotal is not data based, but minLow is)
if doMinPFD        
    % minimize total flux
    % again, to reduce computational burdon, we use the V+ variables in the
    % original MILProblem instead of creating new variables
    %create a new objective function
    % Creating c (objective function)
    c_minFluxLowRxns = ones(length(worm.rxns),1);
    c = [c_ori;c_minFluxLowRxns];
    MILProblem.c = c;
    MILProblem.osense = 1;
    fprintf('Minimizing total flux... \n');
    tic()
    solution = solveCobraMILP_XL(MILProblem, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 0);
    minTotal = solution.obj;
    if solution.stat ~= 1
        error('infeasible or violation occured!');
    end
else
    minTotal = 0;
end
fprintf('Primary Flux Fitting completed! \n');
toc()
% make output of primary flux distribution
FluxDistribution = solution.full(1:length(worm.rxns));
PFD = solution.full(1:length(worm.rxns));
yL = boolean(solution.int((length(Hreactions)+1):(length(Hreactions)+length(Lreactions))));
ClosedLReaction = Lreactions(yL);
OpenedLReaction = Hreactions(boolean(solution.int(1:length(Hreactions))+solution.int(length(Hreactions)+length(Lreactions)+1:2*length(Hreactions)+length(Lreactions))));
N_highFit = length(OpenedLReaction);
N_zeroFit = length(ClosedLReaction);
%% examine the flux distribution
fprintf('the bacterial uptake is: %f\n',FluxDistribution(strcmp(worm.rxns,'EXC0050_L'))); %bac
%fprintf('the ATPm is: %f\n',FluxDistribution(strcmp(worm.rxns,'RCC0005'))); %bac
bacWaste = ...
    FluxDistribution(strcmp(worm.rxns,{'EXC0051_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'protein_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0052_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'rna_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0053_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'dna_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0054_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'peptido_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0055_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'lps_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0056_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'lipa_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0057_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'outerlps_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0058_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'pe_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0059_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'pg_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0060_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'ps_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0063_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'clpn_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0065_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'glycogen_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0072_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'soluble_BAC_L[e]')})+...
    FluxDistribution(strcmp(worm.rxns,{'EXC0142_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'phospholipid_BAC_L[e]')});
wasteDW = bacWaste / (-FluxDistribution(strcmp(worm.rxns,'EXC0050_L')) * MolMass(model.metFormulas{strcmp(model.mets,'BAC_L[e]')}));
fprintf('the total waste of bulk bacteria is: %f\n',wasteDW);
end
