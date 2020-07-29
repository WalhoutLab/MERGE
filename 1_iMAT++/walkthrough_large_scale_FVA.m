%% This walkthrough will guide user to perform iMAT++ combined with FVA analysis for multiple conditions
% We recommend users to first go through the walkthrough_generic.m before starting this walkthrough 
% This will guide user to generate the OFD, upper and lower bounds, and
% the parsed FVA level table for blocking reactions in advanced FPA, for each queried
% conditions.

%% PART I: THE APPLICATION TO THE GENERIC C. ELEGANS MODEL
% We continue to use the RNA-seq data from Bulcha et al, Cell Rep (2019) as an example
% And we use generic C. elegans model for demo purpose. Users could run this script on 
% other models with minor modifications.
% As a demo, we only compare four conditions: N2_OP50, N2_B12, nhr10_OP50 and nhr10_B12

%% Step 1: make the OFDs
% add paths
addpath ~/cobratoolbox/%your cobra toolbox path
addpath /share/pkg/gurobi/810/linux64/matlab/%the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
% Load model
load('iCEL1314.mat');
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake for integration
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
% Prepare flux thresholds (epsilon values)
load('input/wormGeneric/epsilon_generic.mat'); % see walkthrough_generic.m for guidance on generating the epsilon values
conditions = {'N2_OP50', 'N2_B12', 'nhr10_OP50','nhr10_B12'};
% Set parameters
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
%
% First calculate and save the OFDs of all these conditions
for i = 1:length(conditions)
    sampleName = conditions{i};
    load(['input/wormGeneric/exampleGeneCategories/categ_',sampleName,'.mat']) % load categories
    myCSM = struct(); %myCSM: my Context Specific Model
    [myCSM.OFD,...
    myCSM.PFD,...
    myCSM.N_highFit,...
    myCSM.N_zeroFit,...
    myCSM.minLow,...
    myCSM.minTotal,...
    myCSM.minTotal_OFD,...
    myCSM.MILP,...
    myCSM.MILP_PFD,...
    myCSM.HGenes,...
    myCSM.RLNames,...
    myCSM.OpenGene,...
    myCSM.latentRxn,...
    myCSM.Nfit_latent,...
    myCSM.wasteDW]...
    = IMATplusplus(model,epsilon_f,epsilon_r, ExpCateg, modelType);
    save(['output/genericModelDemo/',sampleName,'.mat'],'myCSM');
    eval([sampleName,' = myCSM;']);
    fprintf('OFD of %s is done!\n',sampleName);
end

%% Step 2: FVA
% The FVA calculation is computationally intensive. We expect the speed of ~100
% reactions per miniute in a 4-core Mac laptop. Therefore, the total FVA
% computation time may be ~2 hour in a laptop. For best practice, we
% recommend running it on a multi-core lab server (>=20).
%
% Define the par cluster
myCluster = parcluster('local');
myCluster.NumWorkers = 128;
saveProfile(myCluster);
parpool(20,'SpmdEnabled',false);% adjust according to your computing environment
for i = 1:length(conditions)
    sampleName = conditions{i};
    eval(['myCSM = ',sampleName,';']);
    % setup FVA inputs
    targetRxns = model.rxns;
    parforFlag = 1;
    myFVA = struct(); % save FVA in a structure variable
    % run FVA by calling:
    [myFVA.lb, myFVA.ub] = FVA_MILP(myCSM.MILP, model, targetRxns,parforFlag);
    eval([sampleName,'_FVA = myFVA;']);
    save(['output/genericModelDemo/FVA/',sampleName,'.mat'],'myFVA');
    fprintf('FVA of %s is done!\n',sampleName);
end

%% Step 3: generate reaction status for FPA (level tables)
% We will convert the FVA boundaries to different levels that indicates the
% active/inactive status. The level 1 means "carry flux in OFD", 0 means 
% "not carry flux in OFD, but in ALT" and -1 means "not carry flux in SLNS". 
% Therefore, the metabolic network for a condition consists of reactions
% with levels 0 and 1. Please refer to the paper for more information. 
for z = 1:length(conditions)
    eval(['myCSM = ',conditions{z},';']);
    eval(['myFVA = ',conditions{z},'_FVA;']);
    levels_f = -1*ones(length(model.rxns),1);
    levels_r = -1*ones(length(model.rxns),1);
    % the following decision making algorithm follows the description in
    % Appendix `Flux Variability Analysis of Tissue-Level Metabolic Network
    % Models`
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
% Consistent with the walkthrough of iMAT++, we provide the example of 
% integrating RNA-seq data of 17 human tissues with human model, RECON2.2 
% (Swainston et al., 2016)

% IMPORTANT NOTICE:
% THE FOLLOWING COMPUTATION TAKES ~10 HOURS IN A 20-CORE LAB SERVER. SO, IF
% YOU ARE RUNNING IT ON A LAPTOP, YOU MAY EXPECT LONGER TIME FOR THE PROGRAM
% TO FINISH!

%% step 1: make the OFDs (COMPUTATION SKIPPED)
% Please refer to walkthrough_generic.m on making OFDs for these 17
% tissues. Here we directly load the outputs

% load the model and other inputs
% Add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
% Load original recon2.2 model
load('./input/humanModel/Recon2_2.mat');
% fix a typo in the model
model.genes(strcmp(model.genes,'HGNC:HGNC:2898')) = {'HGNC:2898'};
model.genes(strcmp(model.genes,'HGNC:HGNC:987')) = {'HGNC:987'};
% Constrain the model
model = defineConstraints_iMATpp(model, 1000,0.005);
% parseGPR takes a significant amount of time, so preparse and integrate with the model
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
model = buildRxnGeneMat(model); % some standard fields are missing in the original model, so we generate them
model = creategrRulesField(model);
model_ori = model;
% load epsilons (will be used later)
load('input/humanModel/epsilon.mat'); 
% Remove the dead reactions (that cannot carry flux)
rxns_ori = model.rxns;
model = removeRxns(model,model.rxns(capacity_f == 0 & capacity_r == 0));
[A B] = ismember(model.rxns,rxns_ori);
epsilon_f = epsilon_f(B(A));
epsilon_r = epsilon_r(B(A));
% load the pre-calculated OFDs for 17 tissues (see walkthrough_generic.m)
load('output/humanTissue/outputCollections_NX.mat');
%% Step 2: calculate FVA
% define the par cluster (skip if already done)
% myCluster = parcluster('local');
% myCluster.NumWorkers = 128;
% saveProfile(myCluster);
% parpool(20,'SpmdEnabled',false);% adjust according to your computing environment

targetRxns = model.rxns; % we calculate the FVA of all reactions
% loop through each tissue to calculate the FVA
for i = 1:length(ExampleTissues)
    fprintf('now starting to calculate for %s... \n',ExampleTissues{i});
    fitTime = tic();
    
    % setup parameters for FVA
    % the parameters are important. They are related to the three speed
    % levels discussed in the paper (appendix). Please also check the 
    % technical notes at the end of this walkthrough for details.
    parforFlag = 1;
    RelMipGap = 1e-3;
    
    myFVA = struct(); % the FVA result is stored in a structure variable
    % the OFD (`myCSM`) used for FVA calculation
    myCSM = outputCollections{strcmp(ExampleTissues,ExampleTissues{i})};
    
    % calculate FVA by calling: 
    [myFVA.lb, myFVA.ub] = FVA_MILP(myCSM.MILP_PFD, model, targetRxns,parforFlag,RelMipGap);
    
    save(['output/humanTissue/FVA/',ExampleTissues{i},'.mat'],'myFVA');
    toc(fitTime);
end

%% Step 3: generate reaction status for FPA (level tables)
% we will convert the FVA boundaries to different levels that indicates the
% active/inactive status. The level 1 means "carry flux in OFD", 0 means 
% "not carry flux in OFD, but in ALT" and -1 means "not carry flux in SLNS".
% Therefore, the metabolic network for a condition consists of reactions
% with levels 0 and 1. Please refer to the paper for more information. 

% SPECIAL TREATMENT FOR RECON2.2
% To speed up iMAT++, we removed reactions that don't carry flux (see
% above), such as dead end reactions. However, these dead-end reactions
% may still be useful in metabolite-centric FPA analysis. For example, the 
% melanin biosynthesis pathway is reconstructed in RECON2.2, but there is 
% no reaction to drain melanin at the end. So, the pathway doesn't carry 
% flux. However, the production potential of melanin could still be 
% analyzed by adding in a demand reaction, as we do in metabolite-centric
% FPA. Therefore, we assign these reactions level "0" to not block them.

for z = 1:length(ExampleTissues)
    myCSM = outputCollections{z};
    load(['output/humanTissue/FVA/',ExampleTissues{z},'.mat'])
    
    % We assign levels for all reactions in the original model, since we use the original model in FPA
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

%% TECHNICAL NOTES FOR FVA
%% speed levels and FVA parameter
% As mentioned in our paper (`robustness and usability of MERGE`), we
% developed three speed levels to achieve the computational efficiency that
% meets user's demand. In iMAT++, the speed level is internally controlled
% by a single parameter, `speedMode`. However, it is not parameterized in 
% FVA because it is simply achieved by entering different inputs.

% The speed level of FVA is tuned by two inputs:
% a. MILPproblem_minFlux
% b. relMipGapTol

% Specifically, we have:
% (1) speed level 1 
%     when speed level 1 (speedMode=1) is used in iMAT++, you may want to
%     keep using speed level 1 for FVA. 
%   ```
%     MILPproblem_minFlux = myCSM.MILP (last-step MILP of OFD optimization,
%     see the `MILP` output in iMATplusplus function)
%     relMipGapTol = 1e-12 (default) 
%   ```
%     this setup is used in the FVA demo for C. elegans generic model above
%
% (2) speed level 2
%     when speed level 2 (speedMode=2) is used in iMAT++, you may want to
%     keep using speed level 2 for FVA. 
%   ```
%     MILPproblem_minFlux = myCSM.MILP_PFD (last-step MILP of PFD 
%     optimization, see the `MILP_PFD` output in iMATplusplus function)
%     relMipGapTol = 1e-3
%   ```
%     this setup is used in the FVA demo for human model above
%
% (3) speed level 3
%     when speed level 3 (speedMode=3) is used in iMAT++, you should keep
%     using speed level 3 for FVA. 
%   ```
%     MILPproblem_minFlux = myCSM.MILP_PFD (last-step MILP of PFD 
%     optimization, see the `MILP_PFD` output in iMATplusplus function)
%     relMipGapTol = 1e-3
%   ```
%     The setup of level-3 FVA is the same as level-2 FVA. However, because
%     the level-3 iMAT++ is mathmatically different from level-2 iMAT++,
%     the `MILP_PFD` is not the same. Indeed, the `MILP_PFD` from level-3
%     iMAT++ increases the FVA speed a lot.

% In summary, if you want to run FVA by speed level 1, please follow the
% codes used in FVA demo of generic C. elegans model. If you want to run
% FVA by speed level 2 or 3, please follow the codes in the human model
% demo instead. 






