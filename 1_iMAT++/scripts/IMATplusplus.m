function [OFD,N_highFit,N_zeroFit,minLow,minTotal,OpenGene,wasteDW,HGenes,RLNames,latentRxn,PFD,Nfit_latent,minTotal_OFD,MILP] = IMATplusplus(model,doLatent,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCateg,doMinPFD,latentCAP,modelType,minLowTol,BigModel,verbose)
% Uses the iMAT++ algorithm (`Yilmaz et al., 2020`) to find the optimal flux distribution that fits into a categorized gene expression data. 
% iMAT++ algorithm performs multi-step fitting to find a flux distribution
% that best agrees with rarely, lowly and highly expressed genes in the
% expression profile. It employs a gene-centric feature that minimizes
% false fitting and conflicts due to ambiguous GPR assignment. 
%
% USAGE:
%
%    OFD = IMATplusplus(model,doLatent,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCatag,doMinPFD,latentCAP,modelType)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    doLatent:          whether to perform recursive fitting of latent
%                       reactions
%    storeProp:         (dual-C. elegans tissue network only, use empty ('[]') for other models) the maximum allowed storage molecule uptake, as w/w percentage of bacteria uptake                  
%    SideProp:          (dual-C. elegans tissue network only, use empty ('[]') for other models) the maximum allowed side metabolites uptake, as w/w percentage of bacteria uptake                  
%    epsilon_f:         the epsilon sequence for the forward direction of all reactions. The non-applicable reaction (i.e., for a irreversible,
%                       reverse-only reaction) should have a non-zero default to avoid numeric error.
%    epsilon_r:         the epsilon sequence for the reverse direction of all reactions. The non-applicable reaction (i.e., for a irreversible,
%                       forward-only reaction) should have a non-zero default to avoid numeric error.
%    ATPm:              (C. elegans network only, use empty ('[]') for other models) the Non Growth Associated Maintenance (NGAM) used in fitting               
%    ExpCateg:          Expression Categories used for iMAT++ fitting. The
%                       ExpCateg should be a structure variable with "high", "low", "dynamic" and
%                       "zero" four fields. See the walkthrough scripts on how to generate
%                       them from raw expression quantification data (i.e., TPM)
% OPTIONAL INPUTS:
%    doMinPFD:          whether to perform flux minimization of PFD (after
%                       minimizing the total flux for lowly and rarely expressed genes). PFD
%                       flux minimization is required to obtain OFD, but could be omitted if
%                       one only wants to get MILP fitting result.
%    latentCAP:         the total flux cap for recursive fitting of latent
%                       reactions. The total flux will be capped at (1 +
%                       latentCAP)*OriginalTotalFlux; The default cap is
%                       0.01
%    modelType:         an integer to specify the input COBRA model. 1 ==
%                       dual C. elegans tissue model, 2 == generic C. elegans model, 3 ==
%                       other COBRA models. Optional parameters to control model constraints
%                       like ATP maintenance, side metabolite percentage and etc will be
%                       avaible for setting modelType to 1 and 2. 
%    minLowTol:         the numerical tolerance for total flux of the reactions dependant on lowly and
%                       rarely expressed genes. This parameter tunes the
%                       strigency of constrianing such fluxes to the
%                       minimal level. Default value (1e-5) provides very
%                       strigent constriant that pushes the total flux to
%                       the minimum level. However, in some cases where the
%                       lowly expressed and rarely expressed genes
%                       extensively conflict with highly expressed genes,
%                       we recommand to use a larger tolerance such as the
%                       default epsilon (ideally allowing one mis-fitting)
%    BigModel:          (0 or 1) to indicate if the "big model" mode is
%                       used. This mode is recommanded for all complex models. 
%                       In this mode, the low reactions (dependent on rarely and lowly
%                       expressed genes) will be constrained by rigid boundaries after flux
%                       minimization. The original integer variables and minLow total flux
%                       constriants will be removed. In general, this mode gives almost
%                       identical flux prediction as the normal mode.
%    verbose:           (0 or 1) to show the MILP log or not
%
%
% OUTPUT:
%   OFD:                the Optimal Flux Distribution (OFD)
% OPTIONAL OUTPUTS:
%   N_highFit:          the number of highly expressed genes fitted
%   N_zeroFit:          the number of rarely expressed genes fitted
%   minLow:             the minimal total flux of reactions dependent on
%                       rarely and lowly expressed genes
%   minTotal:           the minimal total flux of PFD (primary flux distribution)
%   OpenGene:           the list of fitted (carrying flux) highly expressed genes
%   wasteDW:            (C. elegans network only) the percentage of wasted bacterial biomass (DW/DW)
%   HGenes:             the list of all the highly expressed genes to be fitted
%   RLNames:            the list of all the reactions dependent on rarely expressed genes
%   latentRxn:          the list of all identified latent reactions to be fitted 
%   PFD:                the primary flux distribution (PFD)
%   Nfit_latent:        the (total) number of latent reactions fitted
%   minTotal_OFD:       the minimal total flux of OFD
%   MILP:               the MILP problem (in COBRA format) in the final
%                       flux minimization of OFD. This MILP can serve as a
%                       startpoint for any custom analysis such as flux minimization or FVA.
%
% `Yilmaz et al. (2020). Final Tittle and journal.
%
% .. Author: - (COBRA implementation) Xuhang Li, Mar 2020
%% this is the main integration function
% apply default parameters
if (nargin < 9)
    doMinPFD = true;
end
if (nargin < 10) % the flux tolerence cap of total flux in the latent fitting. By default, we use 5%
    latentCAP = 0.05;
end
if (nargin < 11) % type of the input model 
    modelType = 1; %1==>dual C.elegans; 2==>generic C. elegans; 3==>non-c.elegans
end
if (nargin < 12) 
    minLowTol = 1e-5; % default value of low flux tolerance 
end
if (nargin < 13) 
    BigModel = 0; % by default, do normal IMAT++
end
if (nargin < 14) 
    verbose = 0; % by default, don't show the MILP details
end
%set global constant 
bacMW=966.28583751; %only will be used for C. elegans model
%changeCobraSolverParams('LP','optTol', 10e-9);
%changeCobraSolverParams('LP','feasTol', 10e-9);
%% mapping the gene categories to reactions
fprintf('Start flux fitting... \n');
tic()
fprintf('Processing expression data... \n');
worm = model;
% process gene expression data
expression_gene=struct;
expression_gene.value = [3*ones(length(ExpCateg.high),1);2*ones(length(ExpCateg.dynamic),1);1*ones(length(ExpCateg.low),1);zeros(length(ExpCateg.zero),1)];
expression_gene.gene = [ExpCateg.high;ExpCateg.dynamic;ExpCateg.low;ExpCateg.zero];
[expressionRxns] = mapExpressionToReactions_xl(worm, expression_gene);
RHNames = worm.rxns(expressionRxns == 3);
RLNames = worm.rxns(expressionRxns == 0); % for the first step integration
% only genes associated with high reactions are high genes!

HGenes = {};
for i = 1: length(ExpCateg.high)
    myrxns = worm.rxns(any(worm.rxnGeneMat(:,strcmp(ExpCateg.high(i),worm.genes)),2),1);
    if any(ismember(myrxns,RHNames))
        HGenes = [HGenes;ExpCateg.high(i)];
    end
end

%% begin to do iMAT++ integration
%% Step1: do MILP integration of highly expressed genes and rarely expressed genes
toc()
fprintf('fitting high genes with MILP... \n');
tic()
% setup model specifc constraints for C. elegans mdoel
if modelType == 1
    if storeProp>=0 && SideProp>=0 && ATPm>=0
        % adjust the ATPm according to the bacterial uptake. This requires the
        % manual tuning
        worm = changeRxnBounds(worm,'RCC0005_I',ATPm,'l');
        worm = changeRxnBounds(worm,'RCC0005_X',ATPm,'l');
        % change the side proportion
        worm.S(end-1, strcmp('EXC0050_L',worm.rxns)) = storeProp*bacMW*0.01;
        worm.S(end, strcmp('EXC0050_L',worm.rxns)) = SideProp*bacMW*0.01;
    else
        error('Please check your input parameter for side proportion, storage proportion, and ATPm!');
    end
elseif modelType == 2
    worm = changeRxnBounds(worm,'RCC0005',ATPm,'l');
elseif modelType == 3
    fprintf('Doing integration for a non-C. elegans model. Make sure you constrained the ATPm before run the integration! \n');
end
% perform the gene-centric MILP fitting of highly expressed genes and reactions dependent on rarely expressed genes
[~, ~,~,solution, MILProblem] = ExpressionIntegrationByMILP(worm, RHNames, RLNames,HGenes, epsilon_f, epsilon_r,[],[],verbose);
toc()
%% Step2: do minimization of low reactions
fprintf('Minimizing low flux... \n');
tic()
% convert the objective value in a solution to the model constraints
MILProblem =  solution2constraint(MILProblem,solution);
NfitInd = length(MILProblem.b);%save the index of this constraint for later use
% add new variables (absolute flux proxy) for flux minimization
MILProblem = addAbsFluxVariables(MILProblem, worm);
% set the initial solution (as the last step solution) to boost speed 
% MILProblem.x0 = [solution.full;abs(solution.full(1:length(worm.rxns)))];
% minimize low flux (reactions dependent on lowly and rarely expressed
% genes)
% please note that, to use minimal variables in the MILP problem, 
% we do the flux minimization by the absolute flux proxy variables we just
% generated for total flux minimization

% create a new objective function
% Creating c vector (objective function)
c_v = zeros(size(worm.S,2),1);
c_yh1 = zeros(length(RHNames),1);
c_yl = zeros(length(RLNames),1);
c_yh2 = zeros(length(RHNames),1);
c_minFluxLowRxns = zeros(length(worm.rxns),1);
c_minFluxLowRxns(expressionRxns == 1 | expressionRxns == 0) = 1; %set the coeffi of "low" and "rare" reactions as 1
c_gene = zeros(length(HGenes),1);
c = [c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns];
MILProblem.c = c;
% set the objective sense as minimize 
MILProblem.osense = 1;
% sovle the MILP problem
solution = solveCobraMILP_XL(MILProblem, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', verbose);%we use customized solver interface to fine-tune the solver parameters for gurobi.
if solution.stat ~= 1
    error('infeasible or violation occured!');
end
toc()
fprintf('Minimizing low flux completed! \n');
minLow = solution.obj;
%% step3: minimize total flux as needed
if BigModel
    % in big model mode, we also fix all the effectively eliminated flux
    MILProblem.ub(expressionRxns == 1 | expressionRxns == 0) = abs(solution.full(expressionRxns == 1 | expressionRxns == 0));
    MILProblem.lb(expressionRxns == 1 | expressionRxns == 0) = -abs(solution.full(expressionRxns == 1 | expressionRxns == 0));
    % then remove redundant variables (the zero reaction integer variables)
    c_v = zeros(size(worm.S,2),1);
    c_yh1 = zeros(length(RHNames),1);
    c_yl = ones(length(RLNames),1);
    c_yh2 = zeros(length(RHNames),1);
    c_minFluxLowRxns = zeros(length(worm.rxns),1);
    c_gene = zeros(length(HGenes),1);
    lowInd = boolean([c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns]);
    MILProblem.A(:,lowInd) = [];%remove variables
    NfitZero = sum(abs(solution.full(expressionRxns == 0)) <= 1e-9);
    MILProblem.b(NfitInd) = MILProblem.b(NfitInd) - NfitZero; %update the Nfit constraints
    MILProblem.lb(lowInd) = [];%update boundary
    MILProblem.ub(lowInd) = [];%update boundary 
    MILProblem.vartype(lowInd) = [];%update variable type
    RLNames_ori = RLNames;
    RLNames = [];%update RLNames (should be empty because variables are deleted)
else
    % set a minFlux tolerance 
    tol = minLowTol; %it could be solver tolarence or some larger number to increase the numeric stability and allow flexible fitting of zero and low reactions
    solution.obj = solution.obj + tol;
    MILProblem = solution2constraint(MILProblem,solution);
end
% MILProblem.x0 = solution.full;
MILP = MILProblem; %write the MILP output. If latent is done, this MILP will be updated.
if doMinPFD
    % minimize total flux
    % again, to reduce computational burdon, we use the same proxy
    % variables used in the low flux minimization
    % create a new objective function
    % Creating c vector (objective function)
    c_v = zeros(size(worm.S,2),1);
    c_yh1 = zeros(length(RHNames),1);
    c_yl = zeros(length(RLNames),1);
    c_yh2 = zeros(length(RHNames),1);
    c_minFluxLowRxns = ones(length(worm.rxns),1);
    c_gene = zeros(length(HGenes),1);
    c = [c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns];
    MILProblem.c = c;
    fprintf('Minimizing total flux... \n');
    tic()
    solution = solveCobraMILP_XL(MILProblem, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', verbose);
    minTotal = solution.obj;
    if solution.stat ~= 1
        error('infeasible or violation occured!');
    end
    solution.obj = solution.obj*(1+latentCAP); % 1% minFlux constraint!
    MILProblem2 = solution2constraint(MILProblem,solution);
    MILProblem2.x0 = solution.full;
    MILP2 = MILProblem2; 
    % this is a special error-handling section. In rare cases, the highly
    % expressed genes are strongly conflicting with rarely expressed genes.
    % so, our minimization of low reaction fluxes will essentially exlude
    % the fitting of most or all high genes. Finally the minimal total flux
    % could be a very small number, even zero. We raise an error for
    % invalid fittings.
    if solution.obj <= minLowTol %minLowTol represents a small flux tolerance 
        error('All the highly expressed genes conflict with at least one rarely expressed gene. The flux fitting is not valid!')
    end
else
    minTotal = 0;
end
fprintf('Primary Flux Fitting completed! \n');
toc()
% make output of primary flux distribution and some optional outputs
FluxDistribution = solution.full(1:length(worm.rxns));
PFD = solution.full(1:length(worm.rxns));
if BigModel %the low fitting is not controled by binary variable
    nameL = worm.rxns(ismember(worm.rxns,RLNames_ori));
    fluxL = solution.full(ismember(worm.rxns,RLNames_ori));
    ClosedLReaction = nameL(abs(fluxL)<1e-9);
else
    nameL = worm.rxns(ismember(worm.rxns,RLNames));
    yL = boolean(solution.int((length(RHNames)+1):(length(RHNames)+length(RLNames))));
    ClosedLReaction = nameL(yL);
end
OpenGene = HGenes(boolean(solution.int((2*length(RHNames)+length(RLNames))+1:end)));
N_highFit = length(OpenGene);
N_zeroFit = length(ClosedLReaction);
if doLatent
    %% step4. make the latent rxns fitting
    [FluxDistribution,latentRxn,Nfit_latent,minTotal_OFD,MILP] = fitLatentFluxes(MILP2, worm,PFD, HGenes,epsilon_f,epsilon_r,latentCAP,verbose);
    OFD = FluxDistribution;
else
    OFD = [];
    latentRxn = [];
    Nfit_latent = [];
    minTotal_OFD = [];
end
%% provide some major evaluation of the flux distribution (only for C. elegans model)
if modelType == 1
    fprintf('the bacterial uptake is: %f\n',FluxDistribution(strcmp(worm.rxns,'EXC0050_L'))); %bac
    %fprintf('the ATPm is: %f\n',FluxDistribution(strcmp(worm.rxns,'RCC0005'))); %bac
    bacWaste = ...
        FluxDistribution(strcmp(worm.rxns,{'EXC0051_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'protein_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0052_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'rna_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0053_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'dna_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0054_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'peptido_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0055_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'lps_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0056_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'lipa_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0057_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'outerlps_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0058_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'pe_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0059_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'pg_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0060_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'ps_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0063_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'clpn_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0065_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'glycogen_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0072_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'soluble_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0142_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'phospholipid_BAC_L[e]')});
    wasteDW = bacWaste / (-FluxDistribution(strcmp(worm.rxns,'EXC0050_L')) * MolMass(model.metFormulas{strcmp(model.mets,'BAC_L[e]')}));
    fprintf('the total waste of bulk bacteria is: %f\n',wasteDW);
elseif modelType == 2
    fprintf('the bacterial uptake is: %f\n',FluxDistribution(strcmp(worm.rxns,'EXC0050'))); %bac
    %fprintf('the ATPm is: %f\n',FluxDistribution(strcmp(worm.rxns,'RCC0005'))); %bac
    bacWaste = ...
        FluxDistribution(strcmp(worm.rxns,{'EXC0051'})) * MolMass(model.metFormulas{strcmp(model.mets,'protein_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0052'})) * MolMass(model.metFormulas{strcmp(model.mets,'rna_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0053'})) * MolMass(model.metFormulas{strcmp(model.mets,'dna_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0054'})) * MolMass(model.metFormulas{strcmp(model.mets,'peptido_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0055'})) * MolMass(model.metFormulas{strcmp(model.mets,'lps_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0056'})) * MolMass(model.metFormulas{strcmp(model.mets,'lipa_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0057'})) * MolMass(model.metFormulas{strcmp(model.mets,'outerlps_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0058'})) * MolMass(model.metFormulas{strcmp(model.mets,'pe_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0059'})) * MolMass(model.metFormulas{strcmp(model.mets,'pg_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0060'})) * MolMass(model.metFormulas{strcmp(model.mets,'ps_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0063'})) * MolMass(model.metFormulas{strcmp(model.mets,'clpn_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0065'})) * MolMass(model.metFormulas{strcmp(model.mets,'glycogen_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0072'})) * MolMass(model.metFormulas{strcmp(model.mets,'soluble_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0142'})) * MolMass(model.metFormulas{strcmp(model.mets,'phospholipid_BAC[e]')});
    wasteDW = bacWaste / (-FluxDistribution(strcmp(worm.rxns,'EXC0050')) * MolMass(model.metFormulas{strcmp(model.mets,'BAC[e]')}));
    fprintf('the total waste of bulk bacteria is: %f\n',wasteDW);
else
    wasteDW = [];
    fprintf('Please inspect the flux distribution manually for a non-C. elegans model\n');
end
%% fix the output for big model mode
if BigModel
    RLNames = RLNames_ori;%recover RLNames list for output
end
end
