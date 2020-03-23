%% this walkthrough guidance will take user through the application of IMAT++ on a generic model of C. elegans, or other metabolic models (we used human model recon2.2 as an example).
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
% the loaded model is already constrained with default constraints. One can
% further modify it as followings:
model = changeRxnBounds(model,'EXC0050',-1,'l');% we allow 1 unit of bacteria for epsilon calculation
% parseGPR takes hugh amount of time, so preparse and integrate with model
% here 
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
%% load the gene expression data
% For making the gene category file from raw expression quantification
% (i.e, TPM), please refer to "./scripts/makeGeneCategories.m" (run "open makeGeneCategories"). Here we
% directly load the premade gene categories
load('input/exampleGeneCategories/categ_N2_OP50.mat')
% Please note the variable naminclature difference: the high refers to
% "highly expressed genes" in the paper, "dynamic" to "moderately
% expressed", "low" to "lowly expressed" and "zero" to "rarely expressed"
%% prepare epsilons
% users can supply their own epsilon sequence for their own purpose. The
% epsilons should be supplies in the order of reactions in the model, and
% should be equal length of the reactions. The forward direction and
% reverse direction should be supplied seperately.
% we provide an epsilon generator following the methods described in the
% paper
[epsilon_f, epsilon_r] = makeEpsilonSeq(model, model.rxns, 0.01, 0.5);
% NOTE: this may take ~5 mins
%% run the integration function 
% we reset some constraints to make the model ready for integration 
% release the input constraints for integration 
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

%% PART II: THE APPLICATION TO ANY METABOLIC MODEL
% applying IMAT++ to other models is not much different from above.
% However, attentions need to be paid to inputs preparation to make sure it
% is in the correct format. Here we provide an example of integrating
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
% we need to constrian the human model according to the media composition
% luckily, the MFA measurement of exchange fluxes of NCI-60 is available
% (Jain et al, Science, 2012). So, we use the MFA data to define the uptake
% rates

% we first define the uptake rates for calculating epsilons
MFAdata = readtable('input/humanModel/mappedFluxData.xlsx','sheet','cleanData');
model = defineConstriants(model, 1e5,0.01, MFAdata);

% parseGPR takes hugh amount of time, so preparse and integrate with model
% here 
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
model = buildRxnGeneMat(model);
model = creategrRulesField(model);
%% prepare epsilons
% users can supply their own epsilon sequence for their own purpose. The
% epsilons should be supplies in the order of reactions in the model, and
% should be equal length of the reactions. The forward direction and
% reverse direction should be supplied seperately.
% we provide an epsilon generator following the methods described in the
% paper

[epsilon_f, epsilon_r] = makeEpsilonSeq(model, model.rxns, 1, 0.5);
save('input/humanModel/epsilon.mat','epsilon_f','epsilon_r');

load('input/humanModel/epsilon.mat');

% NOTE: calculating epsilon for human model may take 20 mins, so we just
% load the pre-calculated value. One can uncomment the codes and run it again
%% flux fitting for each cell line
% we perform the fitting for 3 different cancer types as a demo. (integration of these 3 cell
% lines runs also relatively fast)
allCells = MFAdata.Properties.VariableNames(6:end);
allCells = allCells(1:2:120);%merge replicates
outputCollections = {};
for i = [1,7,19]
    sampleName = allCells{i};
    %% load the gene expression data
    % For making the gene category file from raw expression quantification
    % (i.e, TPM), please refer to "./scripts/makeGeneCategories.m";
    % Additionally, we provided the original script we used to categorize NCI_60 data,
    % named "./scripts/makeGeneCategories_NCI_60.m", for users' reference.
    load(['input/humanModel/categ_',sampleName,'.mat']);
    % Please note the variable naminclature difference: the high refers to
    % "highly expressed genes" in the paper, "dynamic" to "moderately
    % expressed", "low" to "lowly expressed" and "zero" to "rarely expressed"
    %% run the integration function 
    % we reset the constraints by MFA data of the modeled cell line to make the model ready for integration 
    myMFAdata = MFAdata(:,[1:5,find(ismember(MFAdata.Properties.VariableNames,{sampleName,[sampleName,'_1']}))]);
    myMFAdata.(sampleName) = mean(myMFAdata{:,end-1:end},2);%take the average
    myMFAdata(:,end) = [];
    modelTmp = defineConstriants(model, 1e5,0.01, myMFAdata);
    % set parameters
    doLatent = 1;
    doMinPFD = 1;
    latentCAP = 0.05;
    modelType = 3; % 2 for generic C. elegans model. The default (if not specified) is 1, for the tissue model
    minLowTol = 1; % we allow flexible minimization of lowly and rarely expressed reactions. The violation to minimum total flux could be as large as one epsilon.
    % set the non-applicable parameters to -1 (which will be ignored)
    storeProp = -1;%the storage and side is not applicable for generic model
    SideProp = -1;%the storage and side is not applicable for generic model
    ATPm = -1;

    myCSM = struct(); %myCSM: my Context Specific Model
    try
        [myCSM.OFD,myCSM.N_highFit,myCSM.N_zeroFit,myCSM.minLow,myCSM.minTotal,myCSM.OpenGene,myCSM.wasteDW,myCSM.HGenes,myCSM.RLNames,myCSM.latentRxn,myCSM.PFD,myCSM.Nfit_latent,myCSM.minTotal_OFD,myCSM.MILP] = ...
            IMATplusplus(modelTmp,doLatent,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCateg,doMinPFD,latentCAP,modelType,minLowTol);
    catch ME
        myCSM.errorMessage = getReport(ME);
    end
    % the result is stored in variable "myCSM"
    % For understanding each output field, please see IMATplusplus.m
    % the OFD flux distribution is myCSM.OFD
    outputCollections = [outputCollections;{myCSM}];
end
%% NOTE
% 1. speed concern
% 2. local minimal 