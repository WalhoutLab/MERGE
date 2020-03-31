%% this walkthrough guidance will take user through the application of FPA on a generic model of C. elegans, or other metabolic models (we used human model recon2.2 as an example).
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
%% Part I: THE APPLICATION TO THE GENERIC C. ELEGANS MODEL
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
% we use the RNA-seq data from Bulcha et al, Cell Rep (2019) as an example
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
%NOTE: If a gene is undetected, we will use default value of 0 in the
%calculation. (If a reaction is only associated with undetected genes, it
%will have default penalty (which is 1) in the FPA calculation.)
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
%% PART II: THE APPLICATION TO ANY METABOLIC MODEL
% applying FPA to other models is similair. Like the guidence for iMAT++, 
% here we provide an example of integrating RNA-seq data of NCI-60 cancer 
% cell lines (Reinhold WC et al., 2019) to human model, RECON2.2 (Swainston et al., 2016)
%% 1. load the model and prepare the model
load('./../1_IMAT++/input/humanModel/Recon2_2.mat');
% users may add their own constraints here (i.e., nutriential input
% constraints)
% we define the media condition by the measured fluxes from Jain et al, Science, 2012 
MFAdata = readtable('./../1_IMAT++/input/humanModel/mappedFluxData.xlsx','sheet','cleanData');
% since the flux scale is controled by flux allowance, only the binary
% state (open or closed) of the nutrient uptake matters. So, we allow the
% free uptake for all avaiable nutrient, according to the MFA data
model = defineConstriants(model, 1000,0, MFAdata);%NOTE: This function is different from the same-name function in IMAT++ folder!

% special treatment for recon2.2
% we disassociate the lysozomal ATPase reaction with its genes, because it
% have thousands genes associated but the reaction itself is not very
% meaningful for FBA. Keeping it will cause huge speed problem in GPR
% parsing
model.rules(strcmp(model.rxns,'ATPasel')) = {''};

% parseGPR takes hugh amount of time, so preparse and integrate with model
% here 
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
model = buildRxnGeneMat(model);
model = creategrRulesField(model);
% The FPA analysis requires to pre-parse the GPR and attached it as a field
% in the model. Otherwise parsing GPR in each FPA calculation wastes a lot
% of time. 
%% 2. generate the distance matrix (skip this if you already have the matrix)
% this section will guide users to generate the input for metabolic
% distance calculator from the matlab model.
writematrix(model.S,'distance_inputs/Smatrix_regular.txt');
writecell(model.rxns,'distance_inputs/reactions_regular.txt');
writecell(model.mets,'distance_inputs/metabolites_regular.txt');
writematrix(model.lb,'distance_inputs/LB_regular.txt');
writematrix(model.ub,'distance_inputs/UB_regular.txt');
byProducts = {'co2';'amp';'nadp';'nadph';'ppi';'o2';'nadh';'nad';'pi';'adp';'coa';'atp';'h2o';'h';'gtp';'gdp';'etfrd';'etfox';'crn';'fad';'fadh2'};
% add compartment label to byproducts
byProducts = model.mets(ismember(cellfun(@(x) regexprep(x,'\[.\]$',''),model.mets, 'UniformOutput',false),byProducts));
writecell(byProducts,'distance_inputs/byproducts_regular.txt');
% then you can use the output files in "the distance_inputs" folder for
% calculating distances. Please follow the Distance calculator section in
% Github
%% 3. load the expression files, distance matrix, and other optional inputs
% load expression files
% expression matrix can be in plain text and in any normalized
% quantification metric like TPM or FPKM.
expTbl = readtable('./../1_IMAT++/input/humanModel/log2FPKM.csv');% this is the log2(FPKM+1)
% we use the FPKM to calculate the relative expression levels
expTbl{:,3:end} = 2.^expTbl{:,3:end}-1;
% preprocess the expression table
% to facilate the future use of the expression of many samples, we
% re-organize it into a structure variable.
% the FPA matrix will be in the same order as the master_expression
master_expression = {};%we call this variable "master_expression"
geneInd = ismember(expTbl.GeneID, model.genes); %get the index of genes in the model
for i = 1:size(expTbl,2)-2 %we have 60 samples in the example matrix 
    expression = struct();
    expression.genes = expTbl.GeneID(geneInd);
    expression.value = expTbl{geneInd,i+2};
    master_expression{i} = expression;
end

% load the distance matrix
% users can uncomment the following codes to load from distance calculator
% output; here we load directly from saved matlab variable because of file
% size restriction of GitHub
distance_raw = readtable('./../MetabolicDistance/Output/distanceMatrix_recon2_2.txt','FileType','text','ReadRowNames',true); %we load from the output of the distance calculator. For usage of distance calculator, please refer to the section in Github
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
% In general, we recommend tp set penalty for all Exchange, Demand,
% and Sink reactions to 0 to not penaltize the external reactions. Users 
% may need to interactively tune their special penalties for best flux
% distribution in the FPA calculation
extRxns = model.rxns(cellfun(@(x) ~isempty(regexp(x,'^(EX_|sink_|DM_)', 'once')),model.rxns));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));

% we dont recomand any specific special distance for generic model; In the
% dual model, the special distance was used to discourage using of side
% metabolites. Since side/storage metabolites are not applicable for
% generic model, we don't use any special distance. 
% manualDist = {};

%% 4. run basic FPA analysis
% we only compare the first 10 cell lines for the sake of time
master_expression = master_expression(1:10);
% setup some basic parameters for FPA
n = 1.5;
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
% we perform FPA analysis for two reactions as an example
targetRxns = {'LDH_L','biomass_reaction'};
%The FPA is designed with parfor loops for better speed, so we first initial the parpool
parpool(4)
[FP,FP_solutions] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty);
%NOTE: If a gene is undetected, we will use default value of 0 in the
%calculation. (If a reaction is only associated with undetected genes, it
%will have default penalty (which is 1) in the FPA calculation.)
%% 5. make relative flux potential
relFP_f = nan(size(FP,1),length(master_expression));%flux potential for forward rxns
relFP_r = nan(size(FP,1),length(master_expression));%flux potential for reverse rxns
for i = 1:size(FP,1)
    for j = 1:length(master_expression)
        relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end
%% 6. compare with measured lactate production
Lact = MFAdata(49,:);
predSampleName = expTbl.Properties.VariableNames(3:12);
for i = 1:length(predSampleName)
    exp(i) = mean([Lact.(predSampleName{i}),Lact.([predSampleName{i},'_1'])]);
end
fit = fitlm(exp,relFP_r(1,:));
plot(fit)

%% NOTES FOR RUNNING FPA
%% 1. inspect the flux distribution for reported FP values
% the flux distribution is reported for the irreversible model, in the
% FP_solutions. 
% to get irreversible model
model_irrev = convertToIrreversible(model); % convert to irreversible model
% to inspect the flux distribution, we provide a simple flux tracker 
mytbl = listRxn(model_irrev,FP_solutions{1,1}{2}.full,'lac_L[c]')

%% 2. about gene names and expression data
% we allow letters, numbers, dot, dash and colon in gene names. Any other 
% special symbol needs to be added in the regexp function of line 41 in eval_gpr.m.
% if a gene appears multiple times in the expression table, we use the
% sumation of all levels in the GPR parsing step. 