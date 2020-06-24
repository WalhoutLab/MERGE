%% this walkthrough guidance will take user through the application of FPA on a generic model of C. elegans, or other metabolic models (we used human model recon2.2 as an example).
%
% Please reference to: 
% `xxxxx. (xxx). xxxx
%
% .. Author: - Xuhang Li, March, 2020
%% CAUTION ABOUT MATLAB VERSION
% This demo is developed and tested in MATLAB R2019a version. An earlier
% version may encounter errors with "readtable" and "writematrix"
% functions.
%% add required path 
% **** make sure you do this every time you run the demo!****
% add path for required functions/inputs
addpath ./input/
addpath ./scripts/
addpath ./../input/
addpath ./../bins/
initCobraToolbox(false);
%% Part I: THE APPLICATION TO THE GENERIC C. ELEGANS MODEL
%% 1. prepare the model
load('iCEL1314.mat');

% users may add their own constraints here (i.e., nutritional input constraints)
model = changeRxnBounds(model,'EXC0050',-1000,'l');% we allow unlimited bacteria uptake
model = changeRxnBounds(model,'RCC0005',0,'l');% remove the NGAM 
model.S(ismember(model.mets,{'atp_I[c]','h2o_I[c]','adp_I[c]','h_I[c]','pi_I[c]'}), strcmp('DGR0007_L',model.rxns)) = 0;% remove the energy cost for bacteria digestion

% The FPA analysis requires to pre-parse the GPR and attach it as a field
% in the model. Otherwise parsing GPR in each FPA calculation wastes a lot
% of time. 
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
%% 2. load the expression files, distance matrix, and other optional inputs

% ```load expression files```

% expression matrix can be in plain text and in any normalized
% quantification metric like TPM or FPKM.

% we use the RNA-seq data from Bulcha et al, Cell Rep (2019) as an example
expTbl = readtable('exampleExpression.csv');
% For demo purpose, we only analyze the FPA of four conditions in the
% expression dataset.
conditions = {'N2_OP50', 'N2_B12', 'nhr10_OP50','nhr10_B12'};


% ```preprocess the expression table```

% Since expression tables from different data source may have different
% format, we require user to re-organize your table according to the
% following procedures.

% In FPA, we re-organize the expression table into a variable called
% "master_expression".
% ```REQUIREMENT OF "master_expression"```
% To standardize the format and symbols of different expression matrix, we 
% use a "master_expression" variables as the expression input of FPA. The
% format requirements are as follows:
% (1) it must be a cell arrary of structure variables;
% (2) each structure variable contains the expression information of one
%     condition; and the order of the structure variables in the cell array
%     determines the order of FP values in the output;
% (3) each structure variable must contain two fields: "genes" amd "value".
%     "genes" and "value" should be of equal length.
% (4) Every input "genes" should be measured in all conditions. In other
%     words, it is not recommend to have different "genes" list in
%     different conditions.(although it is allowed by FPA)

% make a new master_expression for these four conditions.
master_expression = {};
% get the index of genes in the model
geneInd = ismember(expTbl.Gene_name, model.genes); 
for i = 1:length(conditions)
    expression = struct();
    expression.genes = expTbl.Gene_name(geneInd);
    expression.value = expTbl.(conditions{i})(geneInd);
    master_expression{i} = expression;
end


% ```load the distance matrix```

% we load from the output of the distance calculator. 
% For usage of distance calculator, please refer to the section in Github
distance_raw = readtable('./../MetabolicDistance/Output/distanceMatrix.txt','FileType','text','ReadRowNames',true); 
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


% ```load the special penalties```

% - We set penalty for all Exchange, Demand, Degradation and Sink reactions
%  to 0, to not penalize the external reactions
manualPenalty = table2cell(readtable('manualPenalty_generic.csv','ReadVariableNames',false,'Delimiter',','));

% other possible inputs:
% we don't recomend any special distance for generic model; In the
% dual model, the special distance was used to discourage the use of side
% metabolites. Since side metabolites are not applicable for generic model,
% we don't use any special distance. 
% manualDist = {};

%% 3. run basic FPA analysis

% The FPA is designed with parfor loops for better speed, so we first 
% initial the parpool
parpool(4)


% ```setup some basic parameters for FPA```

n = 1.5; % distance order
changeCobraSolverParams('LP','optTol', 10e-9); % solver parameter
changeCobraSolverParams('LP','feasTol', 10e-9); % solver parameter


% ```we perform FPA analysis for two reactions as an example```

targetRxns = {'RM04432';'RCC0005'};



% ```Then, run FPA by calling:```

[FP,FP_solutions] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty);


% ```calculate the relative flux potential (rFP)```

% Optionally, we can make relative flux potential (rFP) by the following codes
relFP_f = nan(size(FP,1),length(master_expression));% flux potential for forward rxns
relFP_r = nan(size(FP,1),length(master_expression));% flux potential for reverse rxns
for i = 1:size(FP,1)
    for j = 1:length(master_expression)
        relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end

% NOTICE: If a gene is undetected, this gene is ignored in the GPR parsing
% step, but other detected genes associated with a target reaction will still 
% be used to calculate the penalty.
%% 4. run advanced FPA analysis
% As part of the MERGE package, we recommend user to integrate the result
% of iMAT++ to the FPA analysis. That's saying, to block all reactions that
% don't carry flux in the feasible solution space. These reactions are
% identified by FVA analysis conjoined with IMAT++. Please refer to the
% walkthrough tutorial of IMAT++ and "1_IMAT++/walkthrough_large_scale_FVA.m" 
% for getting the FVA result.

% Assuming the FVA is done, we should have the level table ready for each 
% condition (we calculated the level tables for these four conditions in 
% the FVA walkthrough).


% ```get the block list```

% First, let's merge the level tables for each condition
for i = 1:length(conditions)
    load(['./../1_iMAT++/output/genericModelDemo/FVA/',conditions{i},'levels_.mat']);
    levelTbl_f(:,i) = levels_f;
    levelTbl_r(:,i) = levels_r;
end
% because the FPA is done using the irreversible model, so the rxnID is
% different from the original. We provided a function to get the rxnIDs
% to block according to the level tables.
blockList = getBlockList(model,levelTbl_f,levelTbl_r);


% ```run the FPA again with blocklist input```

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

% we can see that the rFP of nhr10-B12 condition is pushed to 0 when the 
% block list is applied.

% we can also see that the FPA prediction recaptures the repressing of
% Propanoyl-CoA:(acceptor) 2,3-oxidoreductase flux by vitamin B12
% treatment, and the loss of flux activation after nhr10 is deleted.
% (Bulcha et al, Cell Rep (2019), PMID: 30625328) 


%% PART II: THE APPLICATION TO ANY METABOLIC MODEL
% Applying FPA to other models follows the same procedure. Consistent with 
% the guidence for iMAT++, we provide an example of integrating RNA-seq 
% data of human tissue expressions to human model, RECON2.2 (Swainston et 
% al., 2016)

%% initiate the parpool (skip if already initiated)
% The FPA is designed with parfor loops for better speed, so we first 
% initial the parpool
parpool(2); % set according to your computational environment
%% 1. prepare the model
addpath ./input/
addpath ./scripts/
addpath input/
addpath ./../bins/


% ```load model```

% load the original Recon2.2
load('./../1_IMAT++/input/humanModel/Recon2_2.mat');
% fix a typo in the original model
model.genes(strcmp(model.genes,'HGNC:HGNC:2898')) = {'HGNC:2898'};
model.genes(strcmp(model.genes,'HGNC:HGNC:987')) = {'HGNC:987'};


% ```define constraints```

% users may add their own constraints here (i.e., nutritional input)
% we use a set of artificial constraints that allow major nutrients like
% glucose and amino acids. 
model = defineConstriants(model, 1000,0);% NOTE: This function is different from the same-name function in IMAT++ folder!
% in short, we made freely available for all nutrients that are allowed to
% be uptaken, since FPA is independent of quantitative exchange constraints. 
% The flux allowance will limit the flux system. 


% ```special treatments for recon2.2```

% remove parentathsis in the reaction ID (which causes problem in distance calculation)
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');
% we disassociate the lysozomal ATPase reaction with its genes, because it
% has hundreds of genes associated but the reaction itself is not very
% meaningful for FBA/FPA. Keeping it will cause slow speed of GPR parsing
model.rules(strcmp(model.rxns,'ATPasel')) = {''};


% ```pre-parse GPR```

% parseGPR takes hugh amount of time, so preparse and integrate with model
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
model = buildRxnGeneMat(model);% generate other missing fields
model = creategrRulesField(model);
%% 2. (optional) generate the distance matrix (skip this if you already have the matrix)
% this section will guide users to generate the input for metabolic
% distance calculator from a COBRA model.
writematrix(model.S,'distance_inputs/Smatrix_regular.txt');
writecell(model.rxns,'distance_inputs/reactions_regular.txt');
writecell(model.mets,'distance_inputs/metabolites_regular.txt');
writematrix(model.lb,'distance_inputs/LB_regular.txt');
writematrix(model.ub,'distance_inputs/UB_regular.txt');
% user needs to determine their own byproduct metabolites (for details, see Methods) 
byProducts = {'co2';'amp';'nadp';'nadph';'ppi';'o2';'nadh';'nad';'pi';'adp';'coa';'atp';'h2o';'h';'gtp';'gdp';'etfrd';'etfox';'crn';'fad';'fadh2'};% the typical set of byproducts
% add compartment label to byproducts
byProducts = model.mets(ismember(cellfun(@(x) regexprep(x,'\[.\]$',''),model.mets, 'UniformOutput',false),byProducts));
writecell(byProducts,'distance_inputs/byproducts_regular.txt');
% then you can use the output files in folder "the distance_inputs" folder 
% for calculating distances. Please follow the Distance calculator section 
% in Github
%% For an earlier MATLAB version, user may consider the following codes:
% dlmwrite('distance_inputs/Smatrix_regular.txt',full(model.S),',');
% fid=fopen('distance_inputs/reactions_regular.txt','w');
% fprintf(fid,'%s\n',model.rxns{:});
% fclose(fid);
% fid=fopen('distance_inputs/metabolites_regular.txt','w');
% fprintf(fid,'%s\n',model.mets{:});
% fclose(fid);
% dlmwrite('distance_inputs/LB_regular.txt',model.lb,',');
% dlmwrite('distance_inputs/UB_regular.txt',model.ub,',');
% fid=fopen('distance_inputs/byproducts_regular.txt','w');
% fprintf(fid,'%s\n',byProducts{:});
% fclose(fid);
%% 3. load the expression files, distance matrix, and other optional inputs

% ```load expression files```

% expression matrix can be in plain text and in any normalized 
% quantification metric like TPM or FPKM.
expTbl = readtable('./../1_IMAT++/input/humanModel/NX/rna_tissue_consensus_metabolic.csv');% this table provides normalized expression (NX) values for all tissue (metabolic genes only)
% (above expression table is extracted from "RNA consensus tissue gene data" 
% downloaded from www.proteinatlas.org as of Jun-1-2020


% ```preprocess the expression table```

% Since expression tables from different data source may have different
% format, we require user to re-organize your table according to the
% following procedures.

% In FPA, we re-organize the expression table into a cell array of 
% structure variables called "master_expression".
% REQUIREMENT OF "master_expression"
% (1) it must be a cell arrary of structure variables;
% (2) each structure variable contains the expression information of one
%     condition; and the order of the structure variables in the cell array 
%     determines the order of Flux Potential (FP) values in the output;
% (3) each structure variable must contain two fields: "genes" amd "value".
%     "genes" and "value" should be of equal length.
% (4) Each of the input "genes" should be measured in ALL conditions. In 
%     other words, it is not recommended to have different "genes" list in
%     different conditions. 

% makde the master_expression for human tissues
master_expression = {};
conditions = unique(expTbl.Tissue); % conditions to calculate FPA on
% determine the genes that are measured in all conditions (requirement #4)
targetGenes = intersect(expTbl.HGNC(strcmp(expTbl.Tissue,conditions{1})), model.genes); % genes both detected and in the model
for i = 2: length(conditions)
    targetGenes = intersect(expTbl.HGNC(strcmp(expTbl.Tissue,conditions{i})),targetGenes);
end
% make the "master_expression" cell array
for i = 1:length(conditions) 
    expression = struct();
    expression.genes = targetGenes; 
    tissueInd = strcmp(expTbl.Tissue,conditions{i});
    geneLabel = expTbl.HGNC(tissueInd);
    NXval = expTbl.NX(tissueInd);
    % reorder
    [A B] = ismember(targetGenes,geneLabel);
    expression.value = NXval(B(A));
    master_expression{i} = expression;
end


% ```load the distance matrix```

% users can uncomment the following codes to load from the distance 
% calculator output; Here we load directly from saved matlab variable
% because of the file size restriction of GitHub
% if you are loading from the output of distance calculator, use the following line:
% distance_raw = readtable('./../MetabolicDistance/Output/distanceMatrix_recon2_2.txt','FileType','text','ReadRowNames',true); %we load from the output of the distance calculator. For usage of distance calculator, please refer to the section in Github

load('./input/distance_raw_recon2_2.mat');
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


% ```set the special penalties (if desired)```
% *** this could be skipped for quick and simple test of FPA ***

% FPA allows user to avoid the penalty for desired reactions. This could
% be done by setting special penalties (`manualPenalty` parameter). 

% In general, we recommend to avoid the penalty of external reactions in FPA.
% Therefore, we set penalty for all `Exchange`, `Demand`, and `Sink` reactions to 0. 
% Users may need to interactively tune their special penalties for the best flux
% distribution in FPA calculation
extRxns = model.rxns(cellfun(@(x) ~isempty(regexp(x,'^(EX_|sink_|DM_)', 'once')),model.rxns));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));
%% 4. run basic FPA analysis
% we calculate FPA for the same 17 tissues used in iMAT++ demo
% users can calculate their tissues of interest by changing the following 
% list.
ExampleTissues = {'cerebral cortex','spinal cord','midbrain','pons and medulla',...%neuronal tissues (<--> worm neuron/glia)
                'duodenum','colon','small intestine',... % digestive tissues (<--> worm intestine)
                'ovary','seminal vesicle','testis',... % reproductive tissues (<--> worm gonad)
                'heart muscle','smooth muscle','tongue','skeletal muscle',... (<--> worm muscle, pharnyx)
                'skin','kidney','liver'};  %(skin and liver <--> worm hypodermis)

[A B] = ismember(ExampleTissues,conditions);
master_expression = master_expression(B(A));
% the FPA will be reported in the same order as tissues in 'ExampleTissues'

% set up some basic parameters for FPA
n = 1.5; % distance order
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);

%% we perform FPA analysis for four targets as a demo
% (1) isocitrate dehydrogenase reaction (reaction-centric)
% (2) alpha-ketoglutarate dehydrogenase reaction (reaction-centric)
% (3) the cellular demand of melanin (metabolite-centric)
% (4) GABA transporter (metabolite-centric)
%% REACTION-CENTRIC ANALYSIS
%% (1) and (2) 
% reaction-centric analysis is simple, as we don't need to block shortcuts. 
% simply set the target reactions for FPA:
targetRxns = {'ICDHyrm','AKGDm'};
% Finally, run FPA by simply calling:
[FP,FP_solutions] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty);

%% METABOLITE-CENTRIC ANALYSIS
% metabolite-centric analysis can be performed in two ways: 
% (1) on a metabolite sink/demand reaction or 
% (2) on a metabolite transporter reaction.
% we illustrate them individually.
%% (3) the cellular demand of melanin (on a metabolite demand reaction)
% first, let's create the new demand reaction for metabolite-centric analysis of melanin
% (You don't need to add the new reaction if the target demand/sink already exists in the model)
myMet = 'melanin[c]';
myRxnName = ['DMN',myMet];
model_new = addReaction(model,myRxnName,'reactionFormula',[myMet,' -->'],'geneRule', 'NA','printLevel',0);

% then, we need to update the distance matrix because of adding the new demand reaction
[distMat_raw_new,labels_new] = updateDistance(myRxnName,model_new,distMat_raw,labels,max(max(distMat_raw)));   
distMat_new = distMat_raw_new;
for ii = 1:size(distMat_new,1)
    for jj = 1:size(distMat_new,2)
        % take the minimal of the raw distance (see methods for details)
        distMat_new(ii,jj) = min([distMat_raw_new(ii,jj),distMat_raw_new(jj,ii)]);
    end
end

% finally, call FPA function again
targetRxns = {myRxnName};
[FP(3,:),FP_solutions(3,:)] = FPA(model_new,targetRxns,master_expression,distMat_new,labels_new,n, manualPenalty);

% Please note, you may need to block the reactions that could uptake your
% target metabolite (i.e, some exchange reaction). These reactions cause 
% shortcut flux in FPA (see methods for explanation). In this melanin 
% example, we don't need to block since there is no such reaction.

% we will illustrate blocking shortcuts in the next example.
 
%% (4) GABA transporter (on a metabolite transporter reaction)
% we use the canonical sodium-dependent GABA transporter as the target transporter
targetRxns = {'ABUTt4_2_r'}; 

% block shortcuts
shortcutRxns = {'ABUTt2r','GABAVESSEC','GABABGTtc','DM_4abut_n_'}; % other GABA transporters and demand reactions that could interfere with the analysis of the target transporter
model_block = changeRxnBounds(model,shortcutRxns,0,'b');

% to do the analysis of GABA degradation potential, we allow the GABA exchange
model_block = changeRxnBounds(model_block,'EX_4abut_e_',-1000,'l');

% finally, call FPA
[FP(4,:),FP_solutions(4,:)] = FPA(model_block,targetRxns,master_expression,distMat,labels,n, manualPenalty);

% please be advised that you can also analyze the import/export potential
% of other GABA transporters (which we blocked in this case). We only show
% the analysis of one GABA transporter (ABUTt4_2_r) as an example.
%% 5. make relative flux potential (rFP)
relFP_f = nan(size(FP,1),length(master_expression));% flux potential for forward rxns
relFP_r = nan(size(FP,1),length(master_expression));% flux potential for reverse rxns
for i = 1:size(FP,1)
    for j = 1:length(master_expression)
        relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end
% the relative flux potentials (rFP) are stored in two variables, relFP_f and relFP_r. 
% '_f' and '_r' represent forward direction and reverse direction, respectively. 
% Rows of the rFP variables are the four queried reactions; And columns
% are the 17 queried tissues. For example, relFP_f(1,1) is the rFP of the forward
% direction of ICDHym, for the first tissue, cerebral cortex.
%% IMPORTANT NOTICE ON METABOLITE-CENTRIC FPA
% If you are doing large-scale metabolite centric FPA analysis,
% please refer to the corresponding section in "TissueFPA.m" on how to
% programmatically block shortcuts and analyze transporter flux potentials.
% (also on how to avoid mapping expression to penalties repeatedly)
% This demo runs met-centric FPA in a simplified way, so that the shortcuts 
% are blocked manually. 
% Blocking shortcuts programmatically is a model-specific task, as
% different models have different nomenclature system and symbols. 
% Therefore, we can't provide a universal function to do the job. Users may 
% want to engineer our code for worm tissue FPA to do large-scale
% metabolite-centric FPA on their own model.
%% 6. Brief Interpretation 
figure(1)
c = categorical(ExampleTissues);
bar(c,relFP_f(1,:))
title('Isocitrate dehydrogenase (NADP+) - plain FPA')

figure(2)
c = categorical(ExampleTissues);
bar(c,relFP_f(2,:))
title('alpha-ketoglutarate dehydrogenase (NADP+) - plain FPA')

figure(3)
c = categorical(ExampleTissues);
bar(c,relFP_f(3,:))
title('Melanin production - plain FPA')

figure(4)
c = categorical(ExampleTissues);
bar(c,relFP_r(4,:))
title('GABA export - plain FPA')

figure(5)
c = categorical(ExampleTissues);
bar(c,relFP_f(4,:))
title('GABA import - plain FPA')

% INTERPRETATION:
% From the rFP value, we clearly see the superiority of skeletal muscle 
% (and other muscles) in TCA cycle flux potential (via ICDHym and AKGDm),
% and skin in producing melanin. 

% Similarly, the import and export potential of GABA are also in line with 
% our understanding of human physiology. Moreover, FPA is able to distinguish 
% the difference between export and import potential of GABA in different 
% tissues. For example, neuronal tissues are better at exporting (producing)
% GABA than liver and kidney. However, their importing (degrading) potentials
% are similar to that of liver and kidney. This is in line with manual 
% inspection of gene expression data. The GABA synthesis gene, GAD1, is 
% enriched specifically in brain; However, the ABAT gene, which can both  
% synthesize and degrade GABA, is enriched in both liver, kidney and brain. 

% Finally, by inspecting the flux distribution of ICDHyrm FPA (see NOTE#1), 
% we find that FPA successfully integrates the local expression information 
% of TCA cycle (with a shunt) in analyzing the flux potential of isocitrate 
% dehydrogenase. The flux route is pyr --> accoa --> cit --> icit --> akg 
% --> mal --> oaa --> cit. Thus, the FPA successfully made a pathway-level 
% prediction.

%% 7. run advanced FPA analysis
% As part of the MERGE package, we recommend user to integrate the result
% of iMAT++ to the FPA analysis. That's saying, to block all reactions that
% don't carry flux in the feasible solution space. These reactions are
% identified by FVA analysis conjoined with IMAT++. Please refer to 
% walkthrough tutorial of IMAT++ and "1_IMAT++/walkthrough_large_scale_FVA.m" 
% for getting the FVA result.

% Assuming the FVA is done, we should have the level table ready for each 
% condition (we calculated the level tables for these 17 conditions in the 
% FVA walkthrough).


% ```make the to-block reaction list```

% first, let's merge the level tables for each tissue
% tissue name in iMAT++ FVA is modified since space is not allowed in file names
ExampleTissues_name2 = {'cerebralCortex','spinalCord','midbrain','ponsAndMedulla',...%neuronal tissues (<--> worm neuron/glia)
                'duodenum','colon','smallIntestine',... % digestive tissues (<--> worm intestine)
                'ovary','seminalVesicle','testis',... % reproductive tissues (<--> worm gonad)
                'heartMuscle','smoothMuscle','tongue','skeletalMuscle',... (<--> worm muscle, pharnyx)
                'skin','kidney','liver'};  %(skin and liver <--> worm hypodermis)
for i = 1:length(ExampleTissues_name2)
    load(['./../1_iMAT++/output/humanTissue/FVA_levels/',ExampleTissues_name2{i},'_levels.mat']);
    levelTbl_f(:,i) = levels_f;
    levelTbl_r(:,i) = levels_r;
end
% NOTE: tissues in the levelTbl MUST BE in the same order as that in the 
% FPA analysis (i.e, `master_expression`)!

% then, get the list of to-block reactions
% because the FPA is done using the irreversible model, so the rxnID is
% different from the original. We provided a function to get the to-block rxnIDs
blockList = getBlockList(model,levelTbl_f,levelTbl_r);


% ```re-run FPA```

% run the FPA again with the blocklist input
% example (1) and (2)
[FP_adv,FP_solutions_adv] = FPA(model,{'ICDHyrm','AKGDm'},master_expression,distMat,labels,n, manualPenalty,{},max(distMat(~isinf(distMat))),blockList);
% example (3)
[FP_adv(3,:),FP_solutions_adv(3,:)] = FPA(model_new,{'DMNmelanin[c]'},master_expression,distMat_new,labels_new,n, manualPenalty,{},max(distMat(~isinf(distMat))),blockList);
% example (4)
[FP_adv(4,:),FP_solutions_adv(4,:)] = FPA(model_block,{'ABUTt4_2_r'},master_expression,distMat,labels,n, manualPenalty,{},max(distMat(~isinf(distMat))),blockList);
% calculate relative potential
relFP_f_adv = nan(size(FP_adv,1),length(master_expression));%flux potential for forward rxns
relFP_r_adv = nan(size(FP_adv,1),length(master_expression));%flux potential for reverse rxns
for i = 1:size(FP_adv,1)
    for j = 1:length(master_expression)
        relFP_f_adv(i,j) = FP_adv{i,j}(1) ./ FP_adv{i,end}(1);
        relFP_r_adv(i,j) = FP_adv{i,j}(2) ./ FP_adv{i,end}(2);
    end
end

%% 8. interpret the advanced FPA result
figure(6)
c = categorical(ExampleTissues);
bar(c,relFP_f_adv(1,:))
title('Isocitrate dehydrogenase (NADP+) - advanced FPA')

figure(7)
c = categorical(ExampleTissues);
bar(c,relFP_f_adv(2,:))
title('alpha-ketoglutarate dehydrogenase (NADP+) - advanced FPA')

figure(8)
c = categorical(ExampleTissues);
bar(c,relFP_f_adv(3,:))
title('Melanin production - advanced FPA')

figure(9)
c = categorical(ExampleTissues);
bar(c,relFP_r_adv(4,:))
title('GABA export - advanced FPA')

figure(10)
c = categorical(ExampleTissues);
bar(c,relFP_f_adv(4,:))
title('GABA import - advanced FPA')

% INTERPRETATION:
% we can find that the tissue-specific networks blocks the flux potential 
% of melanin production in all tissues except for skin. This provides 
% clearer prediction that melanin is only produced in skin. More 
% importantly, advanced FPA greatly increases the resolution of GABA 
% import/export potential in different tissues. By kicking out tissues 
% whose network doesn't support the production or degredation of GABA, we 
% can distinguish low production/consumption from no production/consumption. 
% And we can exlude tissues who expresses GABA transporter but cannot produce 
% or consume GABA metabolically.

%% TECHNICAL NOTES FOR USING FPA
%% 1. To inspect the flux distribution for reported FP values
% the flux distribution is reported for the irreversible model, and is 
% stored in the "FP_solutions". 

% To inspect the flux distribution, first get irreversible model
model_irrev = convertToIrreversible(model);

% Then, we provide a flux tracker for easy-inspection.
mytbl = listRxn(model_irrev,FP_solutions{1,14}{1}.full,'icit[m]');
% mytbl: column#1: rxnID, col#2: rxn flux, col#3: rxn formula, col#4: flux
% contribution to the queried metabolite.

% you can view the "mytbl" variable for flux going in and out of isocitrate in
% the FPA calculation of skeletal muscle for ICDHym. 

%% 2. notice on gene names and expression data
% (1) Gene Name Rules:
%   (a) we allow letters, numbers, dot, dash and colon in gene names. Any other 
%       symbol needs to be added in the regexp function of line 41 in eval_gpr.m.
%   (b) if a gene name appears multiple times in the expression table, we use the
%       sumation of all expression levels in the GPR parsing step. 

% (2) Expression Data Type:
% By default, we assume that the RNA-seq profile is used as expression
% input. Therefore, we, by default, add a pseudocount of 1 to all the genes,  
% to offset the noise when a gene is lowly expressed. If one would like to 
% use other pseudocount value, or not add pseudocount, please modify
% ./scripts/calculatePenalty.m line 48. Please be advised that we DO NOT
% recommend adding pseudocount for microarray dataset!

% (3) If a gene is undetected, this gene is ignored in the GPR parsing
% step, but other detected genes associated with the same reaction will  
% be still used to calculate the reaction-level penalty.