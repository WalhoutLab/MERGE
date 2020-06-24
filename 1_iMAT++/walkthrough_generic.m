%% This walkthrough will take user through the application of IMAT++ to a generic model of C. elegans, or other metabolic models (we used human model recon2.2 as an example).
%% PART I: THE APPLICATION TO THE GENERIC C. ELEGANS MODEL
%% prepare the model
% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/810/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
% load model
load('iCEL1314.mat');
% this model is already constrained with default constraints. 
% One can further modify it as followings:
model = changeRxnBounds(model,'EXC0050',-1,'l');% we allow 1 unit of bacteria for epsilon calculation
% parseGPR takes hugh amount of time, so preparse and save the result in the model
parsedGPR = GPRparser_xl(model);% Extracting GPR data from the model
model.parsedGPR = parsedGPR;
%% load the gene expression data
% For making the gene category file from raw expression quantification
% (i.e, TPM), please refer to "./scripts/makeGeneCategories.m" (run 
% "open makeGeneCategories"). Here we directly load the premade gene 
% categories

load('input/exampleGeneCategories/categ_N2_OP50.mat')
% Please note the category nomenclature difference: the "high" refers to
% "highly expressed genes" in the paper, "dynamic" to "moderately
% expressed", "low" to "lowly expressed" and "zero" to "rarely expressed"
%% prepare epsilons
% users can supply their own epsilon sequence for their own purpose. The
% epsilons should be supplied in the order of reactions in the model, and
% should be in equal length of that of the reactions. The forward direction 
% and reverse direction should be supplied seperately.
% we provide an epsilon generator following the methods described in the
% paper
[epsilon_f, epsilon_r] = makeEpsilonSeq(model, model.rxns, 0.01, 0.5);
save('input/epsilon_generic.mat','epsilon_f', 'epsilon_r');
% NOTE: this may take ~5 mins
%% run the integration function 
% we reset some constraints to make the model ready for integration 
% release the nutrient constraints for integration 
model = changeRxnBounds(model,'EXC0050',-1000,'l');% free bacteria uptake

% set model type
modelType = 2; % 2 for generic C. elegans model. 
% (The default (if not specified) is 1, for the tissue model)
% set the non-applicable parameters to -1 (which will be ignored)

% for other parameters (i.e, ATPm and latentCAP), we use default values

% run iMAT++ and save outputs in a structure
myCSM = struct(); % myCSM: my Context Specific Model
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
% the result is stored in variable "myCSM"
% For understanding each output field, please see `IMATplusplus.m`
% the OFD flux distribution is myCSM.OFD
%% advanced IMAT++ analysis: Flux Variability Analysis (FVA)
% As introduced in the paper, we further performed FVA analysis to measure 
% the feasible space of each reaction, which in turn provides a set of
% reactions to block in FPA analysis. Users can also perform this analysis 
% for their own dataset/model. However, considering the computational
% intensity, we recommend user to run FVA (of all reactions) on a modern 
% lab server (i.e., >=20 cores, >= 32g mems). For running FVA on all 
% reactions and getting the list of reactions to block, please see 
% walkthrough_large_scale_FVA.m

% Here, we provide a demo for running FVA on a few reactions.
% we can calculate the FVA interval by the `MILP` output in IMAT++
targetRxns = {'BIO0010','BIO0001','BIO0002'};
parforFlag = 0; % whether to run FVA in parallel; we choose "no" for demo
[FVA_lb, FVA_ub] = FVA_MILP(myCSM.MILP, model, targetRxns,parforFlag);

% Together, we show how to calculate the FVA boundaries of three queried 
% reactions for N2_OP50 condition

%% PART II: THE APPLICATION TO ANY METABOLIC MODEL
% applying IMAT++ to other models is not much different from above.
% However, attentions need to be paid to the inputs to make sure they
% are in the correct format. Here we provide an example of integrating
% RNA-seq data of human tissues to human model, RECON2.2 (Swainston et al., 
% 2016)
%% prepare the model
% add paths
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
% load model
load('./input/humanModel/Recon2_2.mat');
% fix a typo in the model
model.genes(strcmp(model.genes,'HGNC:HGNC:2898')) = {'HGNC:2898'};
model.genes(strcmp(model.genes,'HGNC:HGNC:987')) = {'HGNC:987'};

% we need to constrain the human model to calculate epsilons
% As a technical demo, we didn't fine-tune the proper constraint set that
% may represent the human blood environment. We used a set of artificial
% constraints that allows glucose and amino acid as major nutrient sources
model = defineConstriants(model, 1000,0.005);

% parseGPR takes huge amount of time, so preparse and integrate with the 
% model
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
% some standard fields are missing in the original model. We generate them
model = buildRxnGeneMat(model); 
model = creategrRulesField(model);
%% prepare epsilons
% users can supply their own epsilon sequence for their own purpose. The
% epsilons should be supplied in the order of reactions in the model, and
% should be in equal length of the reactions. The forward direction and
% reverse direction should be supplied seperately.
% we provide an epsilon generator following the methods described in the 
% paper

% NOTE: calculating epsilon for human model may take 20 mins (in a laptop), 
% so we just load the pre-calculated value. One can use the following codes 
% to run it again

% [epsilon_f, epsilon_r,capacity_f,capacity_r] = makeEpsilonSeq(model, model.rxns, 0.1, 0.5);
% save('input/humanModel/epsilon.mat','epsilon_f','epsilon_r','capacity_f','capacity_r');

load('input/humanModel/epsilon.mat');
% remove the dead reactions (that cannot carry flux)
rxns_ori = model.rxns;
model = removeRxns(model,model.rxns(capacity_f == 0 & capacity_r == 0));
[A B] = ismember(model.rxns,rxns_ori);
epsilon_f = epsilon_f(B(A));
epsilon_r = epsilon_r(B(A));

%% flux fitting for each tissue
% Please notice that, considering the big size of human model, the integration 
% takes longer than that of the C. elegans model. Each tissue takes around
% ~2 mins in a testing laptop, and few tissues may take longer (i.e, liver
% takes 20 mins).

cateTbl = readtable('input/humanModel/NX/Tcatf_nx_consensus.tsv','FileType','text');
% for the sake of time, we only show the application of MERGE pipeline on
% 17 major tissues.
ExampleTissues = {'cerebralCortex','spinalCord','midbrain','ponsAndMedulla',...% neuronal tissues (<--> worm neuron/glia)
                'duodenum','colon','smallIntestine',... % digestive tissues (<--> worm intestine)
                'ovary','seminalVesicle','testis',... % reproductive tissues (<--> worm gonad)
                'heartMuscle','smoothMuscle','tongue','skeletalMuscle',... % muscle tissues (<--> worm muscle, pharnyx)
                'skin','kidney','liver'};  % skin and other organs (<--> worm hypodermis)
outputCollections = {};
TimeConsumed = [];
for i = 1:length(ExampleTissues)
    sampleName = ExampleTissues{i};
    fprintf('Start to fit %s...\n',sampleName);
    %% load the gene expression data
    % For making the gene category file from raw expression quantification
    % (i.e, TPM), please refer to "./scripts/makeGeneCategories.m" and 
    % the python tools in the repo;
    
    % We premade the categories by the python tool and save results in
    % individual .mat files.
    load(['input/humanModel/NX/categ_',sampleName,'.mat']);
    % Please note the category nomenclature difference: the "high" refers 
    % to "highly expressed genes" in the paper, "dynamic" to "moderately
    % expressed", "low" to "lowly expressed" and "zero" to "rarely expressed"
    %% run the integration function 
    % we reset the constraint of major carbon source to make the model
    % ready for integration (user needs to determine which nutrient(s) they
    % want to release to free)
    modelTmp = model;
    modelTmp.lb(ismember(modelTmp.rxns,{'EX_glc(e)'})) = -1000;% free the main carbon source 
    
    % set parameters
    modelType = 3; % 3 for non-C. elegans model
    
    % The following is the only special parameter for (some) non-C. elegans 
    % model. For few very big models like RECON model, user may tune the  
    % `speedMode` parameter to increase the computational speed. 
    % `speedMode` Level 3 gives highest speed, at the cost of optimization 
    %  stringency; Level 1 is the original iMAT configuration used in the 
    %  paper. Level 2 is somehow in between. We use level 2 in our human 
    % demo. (see iMATplusplus.m for details of different speed modes)
    speedMode = 2;
    
    % we recommend users to first try speed mode 1 for any custom model. 
    % The other two modes are recommended when experiencing slow speed.
    % However, we have shown that this speed tuning doesnt influence
    % vast majority of the discoveries in C. elegans tissue modeling. So,
    % we also recommend users to check out level 2 and 3 speedMode for 
    % their models.
    
    % we use default values for other parameters.
        
    myCSM = struct(); %myCSM: my Context Specific Model
    try
        Tstart = tic;
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
            myCSM.wasteDW,] ...
             = IMATplusplus(modelTmp,epsilon_f,epsilon_r, ExpCateg,modelType,speedMode);
        Tlen = toc(Tstart);
    catch ME
        myCSM.errorMessage = getReport(ME);
        Tlen = NaN;
    end
    % the result is stored in variable "myCSM"
    % For understanding each output field, please see IMATplusplus.m
    % the OFD flux distribution is myCSM.OFD
    TimeConsumed = [TimeConsumed;Tlen];
    outputCollections = [outputCollections;{myCSM}];
end
save(['output/humanTissue/outputCollections_NX.mat'],'outputCollections','ExampleTissues','TimeConsumed');
%% visualize the flux distributions in a heatmap (optional)
fluxdata = [];
for i = 1:length(outputCollections)
    fluxdata = [fluxdata,outputCollections{i}.OFD];
end
% normalize by row (/ max)
fluxdata = fluxdata ./ max(abs(fluxdata),[],2);
fluxdata(isnan(fluxdata)) = 0;
% only keep internal rxns for clustering
model.subSystems = [model.subSystems{:}]';
transportEx = strcmp(model.subSystems,'Transport, extracellular');
ExRxn = findExcRxns(model);
ind = ~transportEx & ~ExRxn;
% make clustergram
distMethod = 'euclidean';
cgo=clustergram(fluxdata(ind,:),'RowLabels',model.rxns(ind),'ColumnLabels',ExampleTissues,'RowPDist',distMethod,'ColumnPDist',distMethod);
c=get(cgo,'ColorMap');
cpr=[c(:,1) c(:,1) c(:,2)];
set(cgo,'ColorMap',cpr);
set(cgo,'Symmetric',true);
set(cgo,'DisplayRange',1);

%% advanced IMAT++ analysis: Flux Variability Analysis (FVA)
% As introduced in the paper, we further performed FVA analysis to measure the
% feasible space of each reaction, which in turn provides a set of
% reactions to block in FPA analysis. Users can also perform this analysis for
% their own dataset/model. However, considering the high computational
% demand, we recommend user to run FVA (of all reactions) on a modern lab server (i.e.,
% >=20 cores, >= 32g mems). For running FVA on all reactions and getting
% the list of reactions to block, please see walkthrough_large_scale_FVA.m

% Here, we provide a demo for running FVA on a few reactions.
% we can calculate the FVA interval by the MILP output in IMAT++

targetRxns = model.rxns([1030,1126]);% some arbitury reactions
parforFlag = 0; % whether to run FVA in parallel; we choose "no" for demo
myCSM = outputCollections{1}; % we pick cerebral Cortex for FVA analysis

% RECOMMENDATIONS ON RUNNING FVA
% the FVA is a computational intensive task. So, it's better to determine
% the FVA configuration according to your demand, computational environment,
% and the model complexity. Here are our recommendations:
% (1) if you are comfortable with speedMode level 1 in iMAT++, do FVA in
%     the same way as our tissue modeling. That is, use 'myCSM.MILP' as your
%     input MILP and use '1e-12' (or your preference) as 'relMipGapTol';
% (2) if you are comfortable with speedMode level 2 or 3 in iMAT++, you 
%     should lower your 'relMipGapTol' to '1e-3'; You will have the choice
%     to run FVA on either 'myCSM.MILP_PFD' or 'myCSM.MILP'. However, we
%     recommend to run FVA on 'myCSM.MILP_PFD' for best speed performance. 
%     (see iMATplusplus.m for details about MILP_PFD and MILP)

% In human application, we use 1e-3 and MILP_PFD.
relMipGapTol = 1e-3;
[FVA_lb, FVA_ub] = FVA_MILP(myCSM.MILP_PFD, model, targetRxns,parforFlag, relMipGapTol);

% Together, we show how to calculate the FVA boundaries of two queried
% reactions for human Cerebral Cortex tissue.


%% TECHNICAL NOTES ON RUNNING IMAT++
%% 1. slow computational speed 
% Although the computational speed of IMAT++ is generally fast, user 
% may experience slow speed in some cases. A common cause is that rarely
% expressed genes (related to low/zero reactions) conflict with
% highly expressed genes, creating a very difficult MILP. It usually
% indicates problems with the input transcriptional profile, where it
% disagrees with the metabolic network reconstruction. We recommend users
% to inspect the conflictions for potential biological discovery. 
% Additionally, complex model may also cause speed problem. See
% the note #3 for recommendations on dealing with complex models.
%% 2. numerical issues 
% Since MILP is a NP-hard problem, even the best solver cannot guarantee to
% find the optimal solution of every MILP. In some cases, user may see solver complaining
% "infeasible model" while the input MILP is clearly feasible. This indicates
% the solver failed to find a feasible solution by its heuristic algorithm. For this, iMAT++ (
% and FVA) can automatically tune the solver parameter when it happens. However, this 
% auto-tune function is only designed for Gurobi. User may need to switch to Gurobi or define their 
% own tuning process for their solvers.
%% 3. running IMAT++ on very complex models (such as RECON series)
% The complex models may experience slow speed in solving the MILP. In addition 
% to running IMAT++ on a high-performance workstation and using "bigModel" mode, 
% we recommend the following model preprocessing before running IMAT++.
% (1) removing inactive metabolic functions (reactions that cannot carry
% flux; we did this above for recon2.2)
% (2) removing reactions dependent on rarely expressed genes (this is
% similar to bigModel mode, but directly removing reactions will more
% effectively increase the speed).
% (3) removing reactions that cannot carry flux when reactions in #2 are removed
% (4) use proper epsilon. Too large or too small epsilon will
% decrease the speed of flux minimization, and even generate unrealistic flux
% distribution. User should evaluate the epsilon choice by the flux burden it generates (i.e.,
% uptake rate of main carbon source).
%
% Normally, the bigModel mode is sufficient to enable MERGE integration of
% dozens of conditions on a large model (i.e., ~5000 rxns). The above
% treatment is only useful for super comprehensive model, such as Recon3D.