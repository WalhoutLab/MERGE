function [OFD,N_highFit,N_zeroFit,minLow,minTotal,OpenGene,wasteDW,HGenes,RLNames,latentRxn,PFD,Nfit_latent,minTotal_OFD,MILP] = autoIntegration_latent(model,doLatent,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCatag,doMinPFD,latentCAP)
%% this is the main integration function
%% annotation is not done yet!
%% inputs
% model ... the model; model needs to be constrained
% storeProp,SideProp ... proportion of sides and storage molecules
% epsilon_f/r ... calculated epsilon vector for every reaction
% ExpCatag ... the expressiom category of each gene; note that the gene
% name should in iCEL format

%apply default
if (nargin < 9)
    doMinPFD = true;
end
if (nargin < 10) %fullsensitivity means sensitivity that includes negative sensitivity score (confident no flux), and due to technical reason, in full sensitivity, -1/+1 score means violation of minTotal flux
    latentCAP = 0.01;
end
%set global constant 
bacMW=966.28583751;
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
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
RHNames = worm.rxns(expressionRxns == 3);
RLNames = worm.rxns(expressionRxns == 0); % for the first step integration
% only genes associated with high reactions are high genes!
HGenes = {};
for i = 1: length(ExpCatag.high)
    myrxns = worm.rxns(any(worm.rxnGeneMat(:,strcmp(ExpCatag.high(i),worm.genes)),2),1);
    if any(ismember(myrxns,RHNames))
        HGenes = [HGenes;ExpCatag.high(i)];
    end
end

%% Step1: do MILP integration
toc()
fprintf('fitting high genes with MILP... \n');
tic()
% adjust the ATPm according to the bacterial uptake. This requires the
% manual tuning
worm = changeRxnBounds(worm,'RCC0005_I',ATPm,'l');
worm = changeRxnBounds(worm,'RCC0005_X',ATPm,'l');
% change the side proportion
worm.S(end-1, strcmp('EXC0050_L',worm.rxns)) = storeProp*bacMW*0.01;
worm.S(end, strcmp('EXC0050_L',worm.rxns)) = SideProp*bacMW*0.01;
[~, ~,~,solution, MILProblem] = ExpressionIntegrationByMILP_NoATP_xl(worm, RHNames, RLNames,HGenes, epsilon_f, epsilon_r);
toc()
% Step2: do minimization of low reactions
%set solution as constraint
fprintf('Minimizing low flux... \n');
tic()
MILProblem = solution2constraint(MILProblem,solution);
MILProblem = addAbsFluxVariables(MILProblem, worm);
%minimize low flux 
%to save computation labor, we do the flux minimization by V+ we generated
%in step 1 matrix 
%create a new objective function
% Creating c (objective function)
c_v = zeros(size(worm.S,2),1);
c_yh1 = zeros(length(RHNames),1);
c_yl = zeros(length(RLNames),1);
c_yh2 = zeros(length(RHNames),1);
c_minFluxLowRxns = zeros(length(worm.rxns),1);
c_minFluxLowRxns(expressionRxns == 1 | expressionRxns == 0) = 1; %set the coeffi of low and zero reactions as 1
c_gene = zeros(length(HGenes),1);
c = [c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns];
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

% step3: minimize total flux as needed
% set a minFlux tolerance 
tol = 1e-5; %it could be 1%
solution.obj = solution.obj + tol;
MILProblem = solution2constraint(MILProblem,solution);
MILP = MILProblem; %write the MILP output, this MILP contains all data-based constraints (minTotal is not data based, but minLow is)
if doMinPFD
    % minimize total flux
    % again, to reduce computational burdon, we use the V+ variables in the
    % original MILProblem instead of creating new variables
    %create a new objective function
    % Creating c (objective function)
    c_v = zeros(size(worm.S,2),1);
    c_yh1 = zeros(length(RHNames),1);
    c_yl = zeros(length(RLNames),1);
    c_yh2 = zeros(length(RHNames),1);
    c_minFluxLowRxns = ones(length(worm.rxns),1);
    c_gene = zeros(length(HGenes),1);
    c = [c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns];
    MILProblem.c = c;
    fprintf('Minimizing total flux... \n');
    tic()
    solution = solveCobraMILP_XL(MILProblem, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 0);
    minTotal = solution.obj;
    if solution.stat ~= 1
        error('infeasible or violation occured!');
    end
    solution.obj = solution.obj*(1+latentCAP); % 1% minFlux constraint!
    MILProblem2 = solution2constraint(MILProblem,solution);
    MILP2 = MILProblem2; 
else
    minTotal = 0;
end
fprintf('Primary Flux Fitting completed! \n');
toc()
% make output of primary flux distribution
FluxDistribution = solution.full(1:length(worm.rxns));
PFD = solution.full(1:length(worm.rxns));
%yH = boolean(solution.int(1:length(RHNames)) + solution.int((length(RHNames)+length(RLNames)+1):(2*length(RHNames)+length(RLNames))));
%nameH = worm.rxns(ismember(worm.rxns,RHNames));
%OpenedHReaction = nameH(yH);
nameL = worm.rxns(ismember(worm.rxns,RLNames));
yL = boolean(solution.int((length(RHNames)+1):(length(RHNames)+length(RLNames))));
ClosedLReaction = nameL(yL);
OpenGene = HGenes(boolean(solution.int((2*length(RHNames)+length(RLNames))+1:end)));
N_highFit = length(OpenGene);
N_zeroFit = length(ClosedLReaction);
if doLatent
    %% perform the sensitivity analysis => discarded for doing safak algorithm!
    %sensitivity = AutoSensitivityAanalysis(worm,ExpCatag,storeProp,SideProp,ATPm,epsilon_f,epsilon_r,capacity_f,capacity_r,true,doFullSensitivity);
    %% make the latent rxns fitting
    [FluxDistribution,latentRxn,Nfit_latent,minTotal_OFD] = fitLatentFluxes(MILP2, worm,PFD, HGenes,epsilon_f,epsilon_r,latentCAP);
    %% filter flux
    OFD = FluxDistribution;
    %OFD = fix(FluxDistribution .* 1e7) ./ 1e7;
else
    OFD = [];
    latentRxn = [];
    Nfit_latent = [];
    minTotal_OFD = [];
end
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
