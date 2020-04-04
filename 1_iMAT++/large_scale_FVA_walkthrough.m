%% This walkthrough will guide user to perform iMAT++ combined with FVA analysis for multiple conditions
%% we recommend users to first go through the walkthrough_generic.m before starting this walkthrough 
% This will guide user to generate the OFD, FVA upper and lower bounds and
% the level table for blocking reactions in advanced FPA, for each queried
% conditions.
% we continue to use the RNA-seq data from Bulcha et al, Cell Rep (2019) as an example
% As a demo, we only compare four conditions: N2_OP50, N2_B12, nhr10_OP50 and nhr10_B12
% the following analysis is performed on generic C. elegans model.
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
parpool(4,'SpmdEnabled',false);%adjust according to your computing environment
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
%% step 3: generate the list of reactions to block in FPA (level table)
% we will convert the FVA boundaries to different levels that indicates the
% active/inactive status. The level 1 means "carry flux in OFD", 0 means "not carry flux in OFD, but in ALT" and
% -1 means "not carry flux in SLNS". Please refer to the paper for more
% information. 
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
        else
            levels_f(i) = -1;%level -1 means "not carry flux in SLNS"
        end

        if -myCSM.OFD(i) > 1e-5
            levels_r(i) = 1;
        elseif -myCSM.OFD(i) > 1e-7 && -myFVA.ub(i) > 1e-7
            levels_r(i) = 1;
        elseif -myFVA.lb(i) > max(epsilon_r(i)-1e-5,1e-5)
            levels_r(i) = 0;
        else
            levels_r(i) = -1;
        end
    end
    save(['output/genericModelDemo/FVA/',conditions{z},'levels_.mat'],'levels_f','levels_r');
    fprintf('level table of %s is saved!\n',conditions{z});
end