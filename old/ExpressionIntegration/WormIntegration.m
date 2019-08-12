%% load model
addpath /Users/xuhangli/Desktop/Walhout_Lab/conjoined_model_project/new_matlab_functions
load('iCEL1311.mat');
worm = addDefaultConstraint(model,'nutritionalFree@1');
worm = changeRxnBounds(worm,'RCC0005',0,'l');
%% laod the gene expression data
fname = 'geneExp.txt';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
GeneExpression = jsondecode(str);
fname = 'geneNameTable.txt';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
GeneNameTable = jsondecode(str);
HighGene = GeneExpression.final_high;
LowGene = GeneExpression.final_low;
ZeroGene = GeneExpression.final_zero;
% convert the name
for i = 1:length(HighGene)
    HighGene(i) = {GeneNameTable.(HighGene{i})};
end
for i = 1:length(LowGene)
    LowGene(i) = {GeneNameTable.(LowGene{i})};
end
for i = 1:length(ZeroGene)
    ZeroGene(i) = {GeneNameTable.(ZeroGene{i})};
end
ModerateGene = worm.genes(~ismember(worm.genes,[HighGene;LowGene;ZeroGene]));

%% process gene expression data
expression_gene=struct;
expression_gene.value = [3*ones(length(HighGene),1);2*ones(length(ModerateGene),1);1*ones(length(LowGene),1);zeros(length(ZeroGene),1)];
expression_gene.gene = [HighGene;ModerateGene;LowGene;ZeroGene];
[expressionRxns parsedGPR] = mapExpressionToReactions_xl(worm, expression_gene);
RHNames = worm.rxns(expressionRxns == 3);
RLNames = worm.rxns(expressionRxns == 0); % for the first step integration
%% Step1: do MILP integration
[epsilon_f,epsilon_r] = makeEpsilonSeq(worm,worm.rxns,0.01,0.5);
alpha = 0.022;
load('ATPlinkedRXNs.mat');
[OpenedGene, OpenedaHReaction,ClosedLReaction,solution, MILProblem] = ExpressionIntegrationByMILP_xl(worm, RHNames, RLNames,HighGene, epsilon_f, epsilon_r, alpha, 'RCC0005',ATPlinkedRXNs);
%% Step2: do minimization of low reactions
%set solution as constraint
MILProblem = solution2constraint(MILProblem,solution);
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
c_minFluxLowRxns(expressionRxns == 1) = 1; %set the coeffi of low reactions as 1
c_gene = zeros(length(HighGene),1);
c = [c_v;c_yh1;c_yl;c_yh2;c_minFluxLowRxns;c_gene];
MILProblem.c = c;
% set the objective sense as minimize 
MILProblem.osense = 1;
% sovle it
solution = solveCobraMILP(MILProblem, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 3);
%% step3: minimize total flux
% set a minFlux tolerance 
tol = 0;
solution.obj = solution.obj + tol;
MILProblem = solution2constraint(MILProblem,solution);
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
c_gene = zeros(length(HighGene),1);
c = [c_v;c_yh1;c_yl;c_yh2;c_minFluxLowRxns;c_gene];
MILProblem.c = c;
solution = solveCobraMILP(MILProblem, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 3);
%% make output
FluxDistribution = solution.full(1:length(worm.rxns));
yH = boolean(solution.int(1:length(RHNames)) + solution.int((length(RHNames)+length(RLNames)+1):(2*length(RHNames)+length(RLNames))));
nameH = worm.rxns(ismember(worm.rxns,RHNames));
OpenedHReaction = nameH(yH);
nameL = worm.rxns(ismember(worm.rxns,RLNames));
yL = boolean(solution.int((length(RHNames)+1):(length(RHNames)+length(RLNames))));
ClosedLReaction = nameL(yL);
OpenGene = HighGene(boolean(solution.int((2*length(RHNames)+length(RLNames))+1:end)));
%% examine the flux distribution
FluxDistribution(strcmp(worm.rxns,'BIO0010')) %BIOMASS
FluxDistribution(strcmp(worm.rxns,'RCC0005')) %ATPm
solution.obj %total flux




%% test any function
load('iCEL1311.mat');
worm = addDefaultConstraint(model,'nutritionalFree@1');
[epsilon_f,epsilon_r] = makeEpsilonSeq(worm,worm.rxns,0.01,0.5);
%%

FBA_linearATPm(worm,0.022, 'RCC0005',ATPlinkedRXNs,'max')
