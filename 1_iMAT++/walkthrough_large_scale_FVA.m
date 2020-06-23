%% This walkthrough will guide user to perform iMAT++ combined with FVA analysis for multiple conditions
%% we recommend users to first go through the walkthrough_generic.m before starting this walkthrough 
% This will guide user to generate the OFD, upper and lower bounds, and
% the parsed FVA level table for blocking reactions in advanced FPA, for each queried
% conditions.
%% PART I: THE APPLICATION TO THE GENERIC C. ELEGANS MODEL
% we continue to use the RNA-seq data from Bulcha et al, Cell Rep (2019) as an example

% And we use generic C. elegans model for demo purpose. Users could run this script on 
% other models with minor modifications.

% As a demo, we only compare four conditions: N2_OP50, N2_B12, nhr10_OP50 and nhr10_B12
%% step 1: make the OFDs
% add paths
addpath ~/cobratoolbox/%your cobra toolbox path
addpath /share/pkg/gurobi/810/linux64/matlab/%the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
% load model
load('iCEL1314.mat');
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
% prepare epsilons
load('input/epsilon_generic.mat'); % see walkthrough_generic.m for guidance on making epsilon
conditions = {'N2_OP50', 'N2_B12', 'nhr10_OP50','nhr10_B12'};
% set parameters
doLatent = 1;
doMinPFD = 1;
latentCAP = 0.05;
ATPm = 10;
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
% set the non-applicable parameters to -1 (which will be ignored)
storeProp = -1;%the storage and side is not applicable for generic model
SideProp = -1;%the storage and side is not applicable for generic model

% first calculate and save the OFDs of all these conditions
for i = 1:length(conditions)
    sampleName = conditions{i};
    load(['input/exampleGeneCategories/categ_',sampleName,'.mat']) % load categories
    myCSM = struct(); %myCSM: my Context Specific Model
    [myCSM.OFD,myCSM.N_highFit,myCSM.N_zeroFit,myCSM.minLow,myCSM.minTotal,myCSM.OpenGene,myCSM.wasteDW,myCSM.HGenes,myCSM.RLNames,myCSM.latentRxn,myCSM.PFD,myCSM.Nfit_latent,myCSM.minTotal_OFD,myCSM.MILP] = ...
        IMATplusplus(model,doLatent,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCateg,doMinPFD,latentCAP,modelType);
    save(['output/genericModelDemo/',sampleName,'.mat'],'myCSM');
    eval([sampleName,' = myCSM;']);
    fprintf('OFD of %s is done!\n',sampleName);
end
%% step 2: calculate FVA
% the FVA calculation is computationally intensive. We expect the speed of ~100
% reactions per miniute in a 4-core Mac laptop. Therefore, the total FVA
% computation time may be ~2 hour in a laptop. For best practice, we recommend to run
% it on a multi-core lab server (>=20).

% define the par cluster
myCluster = parcluster('local');
myCluster.NumWorkers = 128;
saveProfile(myCluster);
parpool(20,'SpmdEnabled',false);% adjust according to your computing environment
for i = 1:length(conditions)
    sampleName = conditions{i};
    eval(['myCSM = ',sampleName,';']);
    targetRxns = model.rxns;
    parforFlag = 1;
    myFVA = struct(); %my context specific model
    [myFVA.lb, myFVA.ub] = FVA_MILP(myCSM.MILP, model, targetRxns,parforFlag);
    eval([sampleName,'_FVA = myFVA;']);
    save(['output/genericModelDemo/FVA/',sampleName,'.mat'],'myFVA');
    fprintf('FVA of %s is done!\n',sampleName);
end
%% step 3: generate the list of reactions to block in FPA (level tables)
% we will convert the FVA boundaries to different levels that indicates the
% active/inactive status. The level 1 means "carry flux in OFD", 0 means "not carry flux in OFD, but in ALT" and
% -1 means "not carry flux in SLNS". Please refer to the paper for more information. 
for z = 1:length(conditions)
    eval(['myCSM = ',conditions{z},';']);
    eval(['myFVA = ',conditions{z},'_FVA;']);
    levels_f = -1*ones(length(model.rxns),1);
    levels_r = -1*ones(length(model.rxns),1);
    for i = 1:length(model.rxns)
        if myCSM.OFD(i) > 1e-5 % 1e-5 is the tol_flux, see supplement text for details
            levels_f(i) = 1;%level 1 means "carry flux in OFD"
        elseif myCSM.OFD(i) > 1e-7 && myFVA.lb(i) > 1e-7 % 1e-7 is the tol_zero, see supplement text for details
            levels_f(i) = 1;
        elseif myFVA.ub(i) > max(epsilon_f(i)-1e-5,1e-5)
            levels_f(i) = 0;%level 0 means "not carry flux in OFD, but in ALT"
        elseif isnan(myFVA.ub(i))
            levels_f(i) = nan;
        else
            levels_f(i) = -1;%level -1 means "not carry flux in SLNS"
        end

        if -myCSM.OFD(i) > 1e-5
            levels_r(i) = 1;
        elseif -myCSM.OFD(i) > 1e-7 && -myFVA.ub(i) > 1e-7
            levels_r(i) = 1;
        elseif -myFVA.lb(i) > max(epsilon_r(i)-1e-5,1e-5)
            levels_r(i) = 0;
        elseif isnan(myFVA.lb(i))
            levels_r(i) = nan;
        else
            levels_r(i) = -1;
        end
    end
    save(['output/genericModelDemo/FVA/',conditions{z},'_levels.mat'],'levels_f','levels_r');
    fprintf('level table of %s is saved!\n',conditions{z});
end

%% PART II: THE APPLICATION TO ANY METABOLIC MODEL
% applying IMAT++ to other models is not much different from above.
% Consistent with the walkthrough of iMAT++, we provide the example of integrating
% RNA-seq data of 17 human tissues to human model, RECON2.2 (Swainston et al., 2016)

%IMPORTANT NOTICE:
%THE FOLLOWING COMPUTATION TAKES ~4-5 HOURS IN A 20-CORE LAB SERVER. SO, IF
%YOU ARE RUNNING IT ON A LAPTOP, YOU MAY EXPECT LONGER TIME FOR THE PROGRAM
%TO FINISH!
%% step 1: make the OFDs (SKIPPED)
% Please refer to walkthrough_generic.m on making OFDs for these 17
% tissues. Here we directly load the outputs

% prepare the model and epsilons
% add paths
addpath ~/cobratoolbox/%your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/%the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
% load model
load('./input/humanModel/Recon2_2.mat');
% fix a typo
model.genes(strcmp(model.genes,'HGNC:HGNC:2898')) = {'HGNC:2898'};
model.genes(strcmp(model.genes,'HGNC:HGNC:987')) = {'HGNC:987'};
% constrain the model
model = defineConstriants(model, 1000,0.005);
% parseGPR takes huge amount of time, so preparse and integrate with the model
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
model = buildRxnGeneMat(model); % some standard fields are missing in the original model. We generate them
model = creategrRulesField(model);
model_ori = model;
% load epsilons
load('input/humanModel/epsilon.mat');
% remove the dead reactions (that cannot carry flux)
rxns_ori = model.rxns;
model = removeRxns(model,model.rxns(capacity_f == 0 & capacity_r == 0));
[A B] = ismember(model.rxns,rxns_ori);
epsilon_f = epsilon_f(B(A));
epsilon_r = epsilon_r(B(A));
%% step 2: calculate FVA
% define the par cluster (skip if already done)
% myCluster = parcluster('local');
% myCluster.NumWorkers = 128;
% saveProfile(myCluster);
% parpool(20,'SpmdEnabled',false);% adjust according to your computing environment

% load the OFDs for 17 tissues
load('output/humanTissue/outputCollections_NX.mat');
targetRxns = model.rxns; % we calculate the FVA of both X tissue and I tissue
for i = 1:length(ExampleTissues)
    fprintf('now starting to calculate for %s... \n',ExampleTissues{i});
    fitTime = tic();
    parforFlag = 1;
    RelMipGap = 1e-3;
    myFVA = struct(); %my context specific model
    myCSM = outputCollections{strcmp(ExampleTissues,ExampleTissues{i})};
    [myFVA.lb, myFVA.ub] = FVA_MILP(myCSM.MILP_PFD, model, targetRxns,parforFlag,RelMipGap);
    save(['output/humanTissue/FVA/',ExampleTissues{i},'.mat'],'myFVA');
    toc(fitTime);
end
%% step 3: generate the list of reactions to block in FPA (level tables)
% we will convert the FVA boundaries to different levels that indicates the
% active/inactive status. The level 1 means "carry flux in OFD", 0 means "not carry flux in OFD, but in ALT" and
% -1 means "not carry flux in SLNS". Please refer to the paper for more information. 

% NOTICE FOR RECON2.2
% To speed up iMAT++, we removed reactions that don't carry flux (see
% aboved), such as dead end reactions. However, these dead-end reactions
% may still be useful in metabolite-centric FPA analysis. So, we assign
% these reactions level "0" to not block them in FPA.

for z = 1:length(ExampleTissues)
    myCSM = outputCollections{z};
    load(['output/humanTissue/FVA/',ExampleTissues{z},'.mat'])
    
    % we assign levels for all reactions in the original model, since we use the original model in FPA
    levels_f = zeros(length(model_ori.rxns),1);
    levels_r = zeros(length(model_ori.rxns),1);
    for i = 1:length(model_ori.rxns) 
        trimedInd = strcmp(model.rxns,model_ori.rxns{i});
        if any(trimedInd)
            if myCSM.OFD(trimedInd) > 1e-5 % 1e-5 is the tol_flux, see supplement text for details
                levels_f(i) = 1;%level 1 means "carry flux in OFD"
            elseif myCSM.OFD(trimedInd) > 1e-7 && myFVA.lb(trimedInd) > 1e-7 % 1e-7 is the tol_zero, see supplement text for details
                levels_f(i) = 1;
            elseif myFVA.ub(trimedInd) > max(epsilon_f(trimedInd)-1e-5,1e-5)
                levels_f(i) = 0;%level 0 means "not carry flux in OFD, but in ALT"
            elseif isnan(myFVA.ub(trimedInd))
                levels_f(i) = nan;
            else
                levels_f(i) = -1;%level -1 means "not carry flux in SLNS"
            end

            if -myCSM.OFD(trimedInd) > 1e-5
                levels_r(i) = 1;
            elseif -myCSM.OFD(trimedInd) > 1e-7 && -myFVA.ub(trimedInd) > 1e-7
                levels_r(i) = 1;
            elseif -myFVA.lb(trimedInd) > max(epsilon_r(trimedInd)-1e-5,1e-5)
                levels_r(i) = 0;
            elseif isnan(myFVA.lb(trimedInd))
                levels_r(i) = nan;
            else
                levels_r(i) = -1;
            end
        end
    end
    save(['output/humanTissue/FVA_levels/',ExampleTissues{z},'_levels.mat'],'levels_f','levels_r');
    fprintf('level table of %s is saved!\n',ExampleTissues{z});
end