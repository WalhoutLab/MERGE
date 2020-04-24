%% This walkthrough guidance will take user through the application of IMAT++ on a generic model of C. elegans, or other metabolic models (we used human model recon2.2 as an example).
%% PART I: THE APPLICATION TO THE GENERIC C. ELEGANS MODEL
%% prepare the model
% add paths
addpath ~/cobratoolbox/%your cobra toolbox path
addpath /share/pkg/gurobi/810/linux64/matlab/%the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
% load model
load('iCEL1314.mat');
% the loaded model is already constrained with default constraints. 
% One can further modify it as followings:
model = changeRxnBounds(model,'EXC0050',-1,'l');% we allow 1 unit of bacteria for epsilon calculation
% parseGPR takes hugh amount of time, so preparse and save result in model
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
%% load the gene expression data
% For making the gene category file from raw expression quantification
% (i.e, TPM), please refer to "./scripts/makeGeneCategories.m" (run "open makeGeneCategories"). Here we
% directly load the premade gene categories
load('input/exampleGeneCategories/categ_N2_OP50.mat')
% Please note the category nomenclature difference: the "high" refers to
% "highly expressed genes" in the paper, "dynamic" to "moderately
% expressed", "low" to "lowly expressed" and "zero" to "rarely expressed"
%% prepare epsilons
% users can supply their own epsilon sequence for their own purpose. The
% epsilons should be supplied in the order of reactions in the model, and
% should be in equal length as that of the reactions. The forward direction and
% reverse direction should be supplied seperately.
% we provide an epsilon generator following the methods described in the
% paper
[epsilon_f, epsilon_r] = makeEpsilonSeq(model, model.rxns, 0.01, 0.5);
save('input/epsilon_generic.mat','epsilon_f', 'epsilon_r');
% NOTE: this may take ~5 mins
%% run the integration function 
% we reset some constraints to make the model ready for integration 
% release the nutrient constraints for integration 
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake
% set parameters
doLatent = 1;
doMinPFD = 1;
latentCAP = 0.05;
ATPm = 10;
modelType = 2; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
% set the non-applicable parameters to -1 (which will be ignored)
storeProp = -1;%the storage and side is not applicable for generic model
SideProp = -1;%the storage and side is not applicable for generic model

myCSM = struct(); %myCSM: my Context Specific Model
[myCSM.OFD,myCSM.N_highFit,myCSM.N_zeroFit,myCSM.minLow,myCSM.minTotal,myCSM.OpenGene,myCSM.wasteDW,myCSM.HGenes,myCSM.RLNames,myCSM.latentRxn,myCSM.PFD,myCSM.Nfit_latent,myCSM.minTotal_OFD,myCSM.MILP] = ...
    IMATplusplus(model,doLatent,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCateg,doMinPFD,latentCAP,modelType);
% the result is stored in variable "myCSM"
% For understanding each output field, please see IMATplusplus.m
% the OFD flux distribution is myCSM.OFD
%% advanced IMAT++ analysis: Flux Variability Analysis (FVA)
% As introduced in the paper, we further performed FVA analysis to measure the
% feasible space of each reaction, which in turn provides a set of
% reactions to block in FPA analysis. Users can also perform this analysis for
% their own dataset/model. However, considering the computational
% intensity, we recommend user to run FVA (of all reactions) on a modern lab server (i.e.,
% >=20 cores, >= 32g mems). For running FVA on all reactions and getting
% the list of reactions to block, please see walkthrough_large_scale_FVA.m

% Here, we provide a demo for running FVA on a few reactions.
% we provided two ways of running FVA. 
% we can calculate the FVA interval by the MILP output in IMAT++
targetRxns = {'BIO0010','BIO0001','BIO0002'};
parforFlag = 0; % whether to run FVA in parallel; we choose "no" for demo
[FVA_lb, FVA_ub] = FVA_MILP(myCSM.MILP, model, targetRxns,parforFlag);

% alternatively, the FVA could be a standalone analysis. Users could supply the same
% input for IMATplusplus to the following function, to get the boudaries of
% queried reactions. This is useful when you only want to know the
% active/inactive status of a set of reactions.
minLowTol = 1e-5;
[myFVA.lb, myFVA.ub] = IMATplusplus_FVA(model,doLatent,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCateg,doMinPFD,latentCAP,modelType,minLowTol,targetRxns,parforFlag);

% Together, we show how to calculate the FVA boundaries of three queried reactions for N2_OP50 condition

%% PART II: THE APPLICATION TO ANY METABOLIC MODEL
% applying IMAT++ to other models is not much different from above.
% However, attentions need to be paid to inputs to make sure they
% are in the correct format. Here we provide an example of integrating
% RNA-seq data of NCI-60 cancer cell lines (Reinhold WC et al., 2019) to
% human model, RECON2.2 (Swainston et al., 2016)
%% prepare the model
% add paths
addpath ~/cobratoolbox/%your cobra toolbox path
addpath /share/pkg/gurobi/810/linux64/matlab/%the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
% load model
load('./input/humanModel/Recon2_2.mat');
% we need to constrain the human model according to the media composition
model = defineConstriants(model, 1000,0.005);

% parseGPR takes huge amount of time, so preparse and integrate with the model
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
model = buildRxnGeneMat(model); % some standard fields are missing in the original model. We generate them
model = creategrRulesField(model);
%% prepare epsilons
% users can supply their own epsilon sequence for their own purpose. The
% epsilons should be supplied in the order of reactions in the model, and
% should be in equal length of the reactions. The forward direction and
% reverse direction should be supplied seperately.
% we provide an epsilon generator following the methods described in the paper

% NOTE: calculating epsilon for human model may take 20 mins (in a laptop), so we just
% load the pre-calculated value. One can use the following codes to run it again

% [epsilon_f, epsilon_r,capacity_f,capacity_r] = makeEpsilonSeq(model, model.rxns, 0.1, 0.5);
% save('input/humanModel/epsilon.mat','epsilon_f','epsilon_r','capacity_f','capacity_r');

load('input/humanModel/epsilon.mat');
% remove the dead reactions (that cannot carry flux)
rxns_ori = model.rxns;
model = removeRxns(model,model.rxns(capacity_f == 0 & capacity_r == 0));
[A B] = ismember(model.rxns,rxns_ori);
epsilon_f = epsilon_f(B(A));
epsilon_r = epsilon_r(B(A));

%% flux fitting for each cell line
% we perform the fitting for 3 random cell lines as a technical demo. 
ExampleCells = {'CO_COLO205','CNS_SF_268','BR_HS578T'};
outputCollections = {};
for i = 1:length(ExampleCells)
    sampleName = ExampleCells{i};
    %% load the gene expression data
    % For making the gene category file from raw expression quantification
    % (i.e, TPM), please refer to "./scripts/makeGeneCategories.m";
    % Additionally, we provided the original script we used to categorize NCI_60 data,
    % named "./scripts/makeGeneCategories_NCI_60.m", for users' reference.
    load(['input/humanModel/categ_',sampleName,'.mat']);
    % Please note the category naminclature difference: the "high" refers to
    % "highly expressed genes" in the paper, "dynamic" to "moderately
    % expressed", "low" to "lowly expressed" and "zero" to "rarely expressed"
    %% run the integration function 
    % we reset the constraint of major carbon source to make the model
    % ready for integration (user needs to determine which nutrient(s) they
    % want to release to free)
    modelTmp = model;
    modelTmp.lb(ismember(modelTmp.rxns,{'EX_glc(e)'})) = -1000;% free the main carbon source 
    
    % set parameters
    doLatent = 1;
    doMinPFD = 1;
    latentCAP = 0.05;
    modelType = 3; % 3 for non-C. elegans model
    
    % The following is the only special parameter for (some) non-C. elegans model.
    % For few very big models like human model, user may use the big model
    % mode in iMAT++ to increase the computational speed. In this mode,
    % the flexible fitting of rarely and lowly expressed genes are changed
    % to rigid fitting (by boudaries)
    bigModel = true;
    % we recommand users to first try normal mode (bigModel = false) for any custom model. The
    % big model mode is recommended when experiencing extreme low
    % speed. This mode uses rigid boundary constriants (ub and lb) instead
    % of flexible objective constraints (total flux and integer Nfit) in eliminating
    % flux in lowly and rarely expressed reactions. But it generally
    % makes near-identical result as the normal mode.
    
    % set the non-applicable parameters to -1 (which will be ignored)
    storeProp = -1;%the storage and side is not applicable for generic model
    SideProp = -1;%the storage and side is not applicable for generic model
    ATPm = -1;
    minLowTol = 0; % not applicable to big model mode; use zero to avoid error
        
    myCSM = struct(); %myCSM: my Context Specific Model
    try
        [myCSM.OFD,myCSM.N_highFit,myCSM.N_zeroFit,myCSM.minLow,myCSM.minTotal,myCSM.OpenGene,myCSM.wasteDW,myCSM.HGenes,myCSM.RLNames,myCSM.latentRxn,myCSM.PFD,myCSM.Nfit_latent,myCSM.minTotal_OFD,myCSM.MILP] = ...
            IMATplusplus(modelTmp,doLatent,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCateg,doMinPFD,latentCAP,modelType,minLowTol,bigModel);
    catch ME
        myCSM.errorMessage = getReport(ME);
    end
    % the result is stored in variable "myCSM"
    % For understanding each output field, please see IMATplusplus.m
    % the OFD flux distribution is myCSM.OFD
    outputCollections = [outputCollections;{myCSM}];
end
%% TECHNICAL NOTES ON RUNNING IMAT++
%% 1. slow computational speed 
% Although the computational speed of IMAT++ is generally fast, user 
% may experience slow speed in some cases. The reason is that rarely
% expressed genes (related to low/zero reactions) conflict with
% highly expressed genes, creating a very difficult MILP. It usually
% indicates problems with the input transcriptional profile, where it
% disagrees with the metabolic network reconstruction. We recommend users
% to inspect the conflictions for potential biological discovery. 
% Additionally, extremely complex model may also cause speed problem. See
% the note #3 for recommendations on dealing with complex models.
%% 2. numerical issues 
% Since MILP is a NP-hard problem, even the best solver cannot guarantee to
% find the optimal solution of every MILP. In some cases, user may see solver complaining
% "infeasible model" while the input MILP is clearly feasible. This indicates
% the solver failed to find a feasible solution by its heuristic algorithm. Users can uncomment
% line 165 and line 220 to enable the pre-defined initial solution, in "IMATplusplus.m".
% But user should aware that the solver may still get stuck in local
% optimum when using this option. We recommend users to avoid running into 
% numerical problems by tuning solvers or trying suggestions provided note #3.
%% 3. running IMAT++ on very complex models (such as RECON3D)
% The complex models may experience slow speed in solving the MILP. In addition 
% to running IMAT++ on a high-performance workstation and using "bigModel" mode, 
% we recommend the following model preprocessing before running IMAT++.
% (1) removing inactive metabolic functions (reactions that cannot carry
% flux; we did this above for recon2.2)
% (2) removing reactions that dependent on rarely expressed genes (this is
% similar to bigModel mode, but directly removing reactions will more
% effectively increase the speed).
% (3) removing reactions that cannot carry flux when reactions in #2 are removed
% (4) use proper (default) epsilon. Too large or too small epsilon will
% decrease the speed of flux minimization, and even generate invalid flux
% distribution. User should evaluate the epsilon choice by the flux burden it generates (i.e.,
% uptake rate of main carbon source).
%
% These preprocessings together generally decrease the computation demand to normal level 
% that can be done in a routine lab server.