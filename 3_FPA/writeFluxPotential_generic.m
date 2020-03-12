% This walkthrough script helps user to learn how to set up the inputs for
% running FPA on a generic iCEL1314 or a given COBRA model (human model
% recon2.2 in this case). 
%
% Please reference to: 
% `xxxxx. (xxx). xxxx
%
% .. Author: - Xuhang Li, March, 2020

% add path for required functions/inputs
addpath ./model
addpath ./programs/2.refinement/
addpath ./programs/3.FluxFitting/
%% run FPA for generic iCEL1314
%% 1. load the model and prepare the model
load('iCEL1314.mat');
% users may add there own constraints here (i.e., nutriential input
% constraints)
model = changeRxnBounds(model,'EXC0050',-1000,'l');%we allow unlimited bacteria uptake
model = changeRxnBounds(model,'RCC0005',0,'l');%remove the NGAM 
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
load('distance_raw.mat'); %we load the pre-made distance matrix directly. For non-c.elegans GEM, please refer to the later section "making distance matrix"
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
% load the special penalties
manualPenalty = table2cell(readtable('manualPenalty_generic.csv','ReadVariableNames',false,'Delimiter',','));
% load the special distance
manualDist = table2cell(readtable('manualDist_generic.csv','ReadVariableNames',false,'Delimiter',','));

%% 3. run flux efficiency
% setup some basic parameters for FPA
maxDist = 31;
n = 1.5;
targetRxns = {'RM04432','RCC0005'};

changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
[fluxEfficiency,fluxEfficiency_plus] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist);

