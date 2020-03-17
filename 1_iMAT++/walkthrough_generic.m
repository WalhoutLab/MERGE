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