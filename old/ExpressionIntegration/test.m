addpath /Users/xuhangli/Desktop/Walhout_Lab/conjoined_model_project/new_matlab_functions
load('iJO1366.mat');
ecoli = iJO1366;
ecoli = changeRxnBounds(ecoli, 'ATPM', 0,'l');
%add sink
[num,txt,all] = xlsread('sinks_to_fix_rxn.xlsx','met_list');
sinksRxns = all(2:end,1);
ecoli = addSinkReactions(ecoli, sinksRxns);
%add constrains; biomass unconstrained(but core is forbid); sinks and non_lb_exchange are 0.1
ecoli = xlsConst2model('./open_exchange.xlsx', ecoli);
ecoli = xlsConst2model('DeletedGenes.xlsx',ecoli);
%generate GPR
ecoli = changeObjective (ecoli, 'BIOMASS_Ec_iJO1366_WT_53p95M'); 
ecoli = generateRules(ecoli);
ecoli = buildRxnGeneMat(ecoli);   
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
%% process gene expression data  --- not is a fake data
[num,txt,all] = xlsread('Yans_gene_in_model.xlsx','geneYans_all');
wormSize_all = num;
GeneNamesFull = all(2:end,1);
[num,txt,all] = xlsread('Yans_gene_in_model.xlsx','geneYans_hit');
wormSize_hit = num;
GeneNameExpressed = all(2:end,1);

expressionMimicforROSproject = ones(size(GeneNameExpressed));
expression_gene=struct;
expression_gene.value = expressionMimicforROSproject;%all set to 1 for Yans Gene in model 
expression_gene.gene = GeneNameExpressed;%all the Yans Gene in model
[expressionRxns parsedGPR] = mapExpressionToReactions_xl(ecoli, expression_gene);
RHNames = ecoli.rxns(expressionRxns == 3);

%test RL(random)
RLNames = {'EX_12ppd__S_e';'EX_15dap_e';'EX_23camp_e';'EX_23ccmp_e';'EX_23cgmp_e';'EX_23cump_e';'EX_23dappa_e'};
%% do MILP integration
epsilon = 10e-6*ones(length(RHNames),1);
alpha = 0.01;
[OpenedGene, OpenedaHReaction,ClosedLReaction,solution] = ExpressionIntegrationByMILP_xl(ecoli, RHNames, RLNames,GeneNameExpressed, epsilon, epsilon, alpha, 'ATPM',ecoli.rxns);

%% test any function
load('iCEL1311.mat');
worm = addDefaultConstraint(model,'nutritionalFree@1');
[epsilon_f,epsilon_r] = makeEpsilonSeq(worm,worm.rxns,0.01,0.5);
%%
load('iCEL1311.mat');
worm = addDefaultConstraint(model,'default');
worm = changeRxnBounds(worm,'RCC0005',0,'l');
[s1 Lproblem] = FBA_linearATPm(worm,0.022, 'RCC0005',ATPlinkedRXNs,'max')
%%
worm = changeRxnBounds(worm,'BIO0010',s1.obj,'l');
worm = changeObjective(worm,'RCC0005');
s2= FBA_linearATPm(worm,0.022, 'RCC0005',ATPlinkedRXNs,'min')
%% test anything 
fname = 'temp.txt';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
GeneExpression = jsondecode(str);

