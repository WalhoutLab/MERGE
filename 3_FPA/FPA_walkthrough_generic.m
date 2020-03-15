% This walkthrough script helps user to learn how to set up the inputs for
% running FPA on a generic iCEL1314 or a given COBRA model (human model
% recon2.2 in this case). 
%
% Please reference to: 
% `xxxxx. (xxx). xxxx
%
% .. Author: - Xuhang Li, March, 2020

% add path for required functions/inputs
addpath ./input/
addpath ./scripts/
addpath ./../input/
addpath ./../bins/
%% run FPA for generic iCEL1314
%% 1. load the model and prepare the model
load('iCEL1314.mat');
% users may add their own constraints here (i.e., nutriential input
% constraints)
model = changeRxnBounds(model,'EXC0050',-1000,'l');%we allow unlimited bacteria uptake
model = changeRxnBounds(model,'RCC0005',0,'l');%remove the NGAM 
model.S(ismember(model.mets,{'atp_I[c]','h2o_I[c]','adp_I[c]','h_I[c]','pi_I[c]'}), strcmp('DGR0007_L',model.rxns)) = 0;%remove the energy cost for bacteria digestion
% The FPA analysis requires to pre-parse the GPR and attached it as a field
% in the model. Otherwise parsing GPR in each FPA calculation wastes a lot
% of time. 
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
%% 2. load the expression files, distance matrix, and other optional inputs
% load expression files
% expression matrix can be in plain text and in any normalized
% quantification metric like TPM or FPKM.
expTbl = readtable('exampleExpression.csv');
% preprocess the expression table
% to facilate the future use of the expression of many samples, we
% re-organize it into a structure variable.
% the FPA matrix will be in the same order as the master_expression
master_expression = {};%we call this variable "master_expression"
geneInd = ismember(expTbl.Gene_name, model.genes); %get the index of genes in the model
for i = 1:size(expTbl,2)-4 %we have 12 samples in the example matrix 
    expression = struct();
    expression.genes = expTbl.Gene_name(geneInd);
    expression.value = expTbl{geneInd,i+4};
    master_expression{i} = expression;
end

% load the distance matrix
distance_raw = readtable('./../MetabolicDistance/Output/distanceMatrix.txt','FileType','text','ReadRowNames',true); %we load from the output of the distance calculator. For usage of distance calculator, please refer to the section in Github
labels = distance_raw.Properties.VariableNames;
labels = cellfun(@(x) [x(1:end-1),'_',x(end)],labels,'UniformOutput',false);
distMat_raw = table2array(distance_raw);
distMat_min = zeros(size(distance_raw,1),size(distance_raw,2));
for i = 1:size(distMat_min,1)
    for j = 1:size(distMat_min,2)
        distMat_min(i,j) = min([distMat_raw(i,j),distMat_raw(j,i)]);
    end
end
distMat = distMat_min;
% load the special penalties - We set penalty for all Exchange, Demand,
% Degradation and Sink reactions to 0 to not penaltize the external
% reactions
manualPenalty = table2cell(readtable('manualPenalty_generic.csv','ReadVariableNames',false,'Delimiter',','));
% we dont recomand any specific special distance for generic model; In the
% dual model, the special distance was used to discourage using of side
% metabolites. Since side/storage metabolites are not applicable for
% generic model, we don't use any special distance. 
% manualDist = {};

%% 3. run basic FPA analysis
% setup some basic parameters for FPA
n = 1.5;
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
% we perform FPA analysis for two reactions as an example
targetRxns = {'RM04432';'RCC0005'};
%The FPA is designed with parfor loops for better speed, so we first initial the parpool
parpool(4)
[fluxEfficiency,fluxEfficiency_plus] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty);
%NOTE: If a gene is undetected, we will assume this reaction to have zero
%expression in the calculation of weight for an associated reaction. 
%% 4. run advanced FPA analysis
% As part of the MERGE package, we recommand user to integrate the result
% of iMAT++ to the FPA analysis. That's saying, to block all reactions that
% don't carry flux in the feasible solution space. (Please refer to "FVA
% analysis of iMAT++" for getting this reaction list)

% assume the FVA is done, we have the list of no-flux reactions 
load('FVAtbl.mat')
% because the FPA is done using the irreversible model, so the rxnID is
% different from the original. We provided a function to get the rxn list
% for blockList input of FPA
blockList = getBlockList(model,FVAtbl);
% run the FPA again
[fluxEfficiency2,fluxEfficiency_plus2] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,{},max(distMat(~isinf(distMat))),blockList);
%% run FPA for human model RECON 2.2
