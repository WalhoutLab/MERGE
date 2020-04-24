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
initCobraToolbox(false);
%% Part I: THE APPLICATION TO THE GENERIC C. ELEGANS MODEL
%% 1. load the model and prepare the model
load('iCEL1314.mat');
% users may add their own constraints here (i.e., nutriential input constraints)
model = changeRxnBounds(model,'EXC0050',-1000,'l');% we allow unlimited bacteria uptake
model = changeRxnBounds(model,'RCC0005',0,'l');% remove the NGAM 
model.S(ismember(model.mets,{'atp_I[c]','h2o_I[c]','adp_I[c]','h_I[c]','pi_I[c]'}), strcmp('DGR0007_L',model.rxns)) = 0;% remove the energy cost for bacteria digestion
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
% NOTE: the FPA matrix will be in the same order as the master_expression!

% For demo purpose, we only analyze the FPA of four conditions in the
% expression dataset.
conditions = {'N2_OP50', 'N2_B12', 'nhr10_OP50','nhr10_B12'};
% make a new master_expression for these four conditions.
master_expression = {};% we call this variable "master_expression"
geneInd = ismember(expTbl.Gene_name, model.genes); % get the index of genes in the model
for i = 1:length(conditions)
    expression = struct();
    expression.genes = expTbl.Gene_name(geneInd);
    expression.value = expTbl.(conditions{i})(geneInd);
    master_expression{i} = expression;
end

% load the distance matrix
distance_raw = readtable('./../MetabolicDistance/Output/distanceMatrix.txt','FileType','text','ReadRowNames',true); % we load from the output of the distance calculator. For usage of distance calculator, please refer to the section in Github
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
% Degradation and Sink reactions to 0, to not penaltize the external reactions
manualPenalty = table2cell(readtable('manualPenalty_generic.csv','ReadVariableNames',false,'Delimiter',','));

% we dont recomand any specific special distance for generic model; In the
% dual model, the special distance was used to discourage using of side
% metabolites. Since side/storage metabolites are not applicable for
% generic model, we don't use any special distance. 
% manualDist = {};

%% 3. run basic FPA analysis
% setup some basic parameters for FPA
n = 1.5; % distance order
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);

% we perform FPA analysis for two reactions as an example
targetRxns = {'RM04432';'RCC0005'};

% The FPA is designed with parfor loops for better speed, so we first initial the parpool
parpool(4)

% Then, run FPA by calling:
[FP,FP_solutions] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty);

% NOTE: If a gene is undetected, this gene is ignored in the GPR parsing
% step, but any detected genes associated with a target reaction will still 
% be used to calculate the penalty.

% optionally, we can make relative flux potential (rFP) by
relFP_f = nan(size(FP,1),length(master_expression));% flux potential for forward rxns
relFP_r = nan(size(FP,1),length(master_expression));% flux potential for reverse rxns
for i = 1:size(FP,1)
    for j = 1:length(master_expression)
        relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end
%% 4. run advanced FPA analysis
% As part of the MERGE package, we recommend user to integrate the result
% of iMAT++ to the FPA analysis. That's saying, to block all reactions that
% don't carry flux in the feasible solution space. These reactions are
% identified by FVA analysis conjoined with IMAT++. Please refer to walkthrough
% tutorial of IMAT++ and "1_IMAT++/walkthrough_large_scale_FVA.m" for getting the FVA result.

% assume the FVA is done, we have the level table for each condition. (we
% calculated the level tables for these four conditions in the FVA walkthrough)

% then, let's merge the level tables for each condition
for i = 1:length(conditions)
    load(['./../1_iMAT++/output/genericModelDemo/FVA/',conditions{i},'levels_.mat']);
    levelTbl_f(:,i) = levels_f;
    levelTbl_r(:,i) = levels_r;
end

% because the FPA is done using the irreversible model, so the rxnID is
% different from the original. We provided a function to get the rxn list
% to block according to the level tables
blockList = getBlockList(model,levelTbl_f,levelTbl_r);

% run the FPA again
[FP_adv,FP_solutions_adv] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,{},max(distMat(~isinf(distMat))),blockList);
% relative potential
relFP_f_adv = nan(size(FP_adv,1),length(master_expression));%flux potential for forward rxns
relFP_r_adv = nan(size(FP_adv,1),length(master_expression));%flux potential for reverse rxns
for i = 1:size(FP_adv,1)
    for j = 1:length(master_expression)
        relFP_f_adv(i,j) = FP_adv{i,j}(1) ./ FP_adv{i,end}(1);
        relFP_r_adv(i,j) = FP_adv{i,j}(2) ./ FP_adv{i,end}(2);
    end
end
%% Compare the advance FPA with basic FPA
figure(1)
c = categorical(regexprep(conditions,'_','-'));
bar(c,relFP_f(1,:))
title('rFP of Propanoyl-CoA:(acceptor) 2,3-oxidoreductase flux')
figure(2)
c = categorical(regexprep(conditions,'_','-'));
bar(c,relFP_f_adv(1,:))
title('advanced rFP of Propanoyl-CoA:(acceptor) 2,3-oxidoreductase flux')
% we can see that the rFP of nhr10-B12 condition is pushed to 0 when the block
% list is applied.




%% PART II: THE APPLICATION TO ANY METABOLIC MODEL
% applying FPA to other models is similair. Consistent with the guidence for iMAT++, 
% here we provide an example of integrating RNA-seq data of NCI-60 cancer 
% cell lines (Reinhold WC et al., 2019) to human model, RECON2.2 (Swainston et al., 2016)
%% 1. prepare the model
load('./../1_IMAT++/input/humanModel/Recon2_2.mat');
% users may add their own constraints here (i.e., nutriential input constraints)
% we define the constraints according to the media composition
model = defineConstriants(model, 1000,0);% NOTE: This function is different from the same-name function in IMAT++ folder!

% remove parentathsis in the reaction ID (which causes problem in distance calculation)
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');

% a special treatment for recon2.2
% we disassociate the lysozomal ATPase reaction with its genes, because it
% has hundreds of genes associated but the reaction itself is not very
% meaningful for FBA. Keeping it will cause huge speed problem in GPR parsing
model.rules(strcmp(model.rxns,'ATPasel')) = {''};

% parseGPR takes hugh amount of time, so preparse and integrate with model
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
model = buildRxnGeneMat(model);% generate other missing fields
model = creategrRulesField(model);
%% 2. generate the distance matrix (skip this if you already have the matrix)
% this section will guide users to generate the input for metabolic
% distance calculator from the COBRA model.
writematrix(model.S,'distance_inputs/Smatrix_regular.txt');
writecell(model.rxns,'distance_inputs/reactions_regular.txt');
writecell(model.mets,'distance_inputs/metabolites_regular.txt');
writematrix(model.lb,'distance_inputs/LB_regular.txt');
writematrix(model.ub,'distance_inputs/UB_regular.txt');
byProducts = {'co2';'amp';'nadp';'nadph';'ppi';'o2';'nadh';'nad';'pi';'adp';'coa';'atp';'h2o';'h';'gtp';'gdp';'etfrd';'etfox';'crn';'fad';'fadh2'};
% add compartment label to byproducts
byProducts = model.mets(ismember(cellfun(@(x) regexprep(x,'\[.\]$',''),model.mets, 'UniformOutput',false),byProducts));
writecell(byProducts,'distance_inputs/byproducts_regular.txt');
% then you can use the output files in folder "the distance_inputs" folder for
% calculating distances. Please follow the Distance calculator section in Github
%% 3. load the expression files, distance matrix, and other optional inputs
% load expression files
% expression matrix can be in plain text and in any normalized quantification metric like TPM or FPKM.
expTbl = readtable('./../1_IMAT++/input/humanModel/log2FPKM.csv');% this is the log2(FPKM+1)
% we use the FPKM to calculate the relative expression levels
expTbl{:,3:end} = 2.^expTbl{:,3:end}-1;

% preprocess the expression table
% To facilate the future use of the expression of many samples, we re-organize it into a structure variable.
% the order of final FP matrix will be in the same order as the master_expression
master_expression = {};% we call this variable "master_expression"
geneInd = ismember(expTbl.GeneID, model.genes); % get the index of genes in the model
for i = 1:size(expTbl,2)-2 % we have 60 samples in the example matrix 
    expression = struct();
    expression.genes = expTbl.GeneID(geneInd);
    expression.value = expTbl{geneInd,i+2};
    master_expression{i} = expression;
end

% load the distance matrix
% users can uncomment the following codes to load from distance calculator
% output; here we load directly from saved matlab variable because of file size restriction of GitHub

% if you are loading from the output of distance calculator, use the following line
% distance_raw = readtable('./../MetabolicDistance/Output/distanceMatrix_recon2_2.txt','FileType','text','ReadRowNames',true); %we load from the output of the distance calculator. For usage of distance calculator, please refer to the section in Github

load('./input/distance_raw_recon2_2.mat');% the original text file is too large for GitHub
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
% In general, we recommend to set penalty for all Exchange, Demand,
% and Sink reactions to 0, to not penaltize the external reactions. Users 
% may need to interactively tune their special penalties for best flux
% distribution in the FPA calculation
extRxns = model.rxns(cellfun(@(x) ~isempty(regexp(x,'^(EX_|sink_|DM_)', 'once')),model.rxns));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));

% we don't recomend any specific special distance for generic model; In the
% dual model, the special distance was used to discourage using of side
% metabolites. Since side/storage metabolites are not applicable for
% generic model, we don't use any special distance. 
% manualDist = {};

%% 4. run basic FPA analysis
% we only compare the first 10 cell lines for the sake of time
master_expression = master_expression(1:10);
% set up some basic parameters for FPA
n = 1.5; % distance order
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);

% we perform FPA analysis for two reactions as an example
targetRxns = {'ICDHyrm','EX_ach_e_'}; 
% we use the isocitrate dehydrogenase reaction (reaction centric) and production of acetylcholine (metabolite centric) as an example

% The FPA is designed with parfor loops for better speed, so we first initial the parpool
parpool(4)

% Finally, run FPA by simply calling:
[FP,FP_solutions] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty);
% NOTE: If a gene is undetected, this gene is ignored in the GPR parsing
% step, but any detected genes associated with a target reaction will still 
% be used to calculate the penalty.
%% 5. make relative flux potential
relFP_f = nan(size(FP,1),length(master_expression));% flux potential for forward rxns
relFP_r = nan(size(FP,1),length(master_expression));% flux potential for reverse rxns
for i = 1:size(FP,1)
    for j = 1:length(master_expression)
        relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end
% the relative flux potentials (rFP) of 'ICDHym' and 'EX_ach_e_' are in
% the relFP_f and relFP_r. Rows are the two queried reactions and columns
% are the 10 queried conditions. For example, relFP_f(1,1) is the rFP of the forward
% direction of ICDHym for the first cell line, BR_MCF7.

%% TECHNICAL NOTES FOR USING FPA
%% 1. inspect the flux distribution for reported FP values
% the flux distribution is reported for the irreversible model, and is in the "FP_solutions". 
% To inspect the flux distribution, first get irreversible model
model_irrev = convertToIrreversible(model); % convert to irreversible model
% the flux is reported in the order of rxns in model_irrev

% addditionally, we provide a simple flux tracker for easy-inspection.
mytbl = listRxn(model_irrev,FP_solutions{1,1}{1}.full,'icit[m]');
% you can view the "mytbl" variable for flux in and out of isocitrate in
% the FPA calculation of BR_MCF7 cell line for ICDHym. 

%% 2. notice on gene names and expression data
% Gene Name Rules:
% we allow letters, numbers, dot, dash and colon in gene names. Any other 
% special symbol needs to be added in the regexp function of line 41 in eval_gpr.m.

% if a gene appears multiple times in the expression table, we use the
% sumation of all levels in the GPR parsing step. 

% Expression Data Notice:
% By default, we assume a regular RNA-seq profile is used as expression
% input. Therefore, we add a pseudocount of 1 to all the genes. This is to 
% offset the noise when a gene is lowly expressed. If one would like to use
% other pseudocount value, or not add pseudocount, please modify
% ./scripts/calculatePenalty.m line 45. Please be advised that we DO NOT
% recommend adding pseudocount for microarray dataset!