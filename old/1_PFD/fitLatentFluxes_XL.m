function [OFD, latentRxn,Nfit_latent,minTotal] = fitLatentFluxes(MILProblem, model, sensitivity, Hgenes, epsilon_f,epsilon_r)
% the input MILP should have Nfit constraint and MinLow constraint already
% (identical to the MILP before final flux minmization)
% tip: the MILP in the final flux minization could be used directly, the
% objective will be changed anyway.
% this is the beta version of latent calculation that use total minFlux
% indenpendent procedure!
%% define candidate reactions from HGene list
% define the global latent candidate that fulfill the expression criteria
latentCandi = {};
for i = 1:length(model.rxns)
    mygenes = model.genes(logical(model.rxnGeneMat(i,:)));
    if ~any(~ismember(mygenes,Hgenes))
        latentCandi = [latentCandi;model.rxns(i)];
    end
end
% define the initial active reaction set
initialActRxns = unique(sensitivity.rxnID(sensitivity.Score == 3 | sensitivity.Score == 2));
fprintf('starting the latent searching loop... \n');
%% start of the latent rxn searching loop
actRxns = initialActRxns;
latentRxn = {};
while 1
    %% define latent reactions
    % mark the active metabolite list
    ActMet = model.mets(any(model.S(:,ismember(model.rxns,actRxns)),2));
    % get the new latent rxns set
    newLatent = {};
    for i = 1:length(latentCandi)
        % new added 0723// considering the reversibility
        metsL = model.mets(model.S(:,strcmp(model.rxns,latentCandi{i}))<0);
        metsR = model.mets(model.S(:,strcmp(model.rxns,latentCandi{i}))>0);
        rxnInd = strcmp(model.rxns,latentCandi{i});
        if ((model.lb(rxnInd)<0 && model.ub(rxnInd)>0) && ((~any(~ismember(metsL,ActMet))) || (~any(~ismember(metsR,ActMet)))))...% reversible & one side mets all have flux
           ||...
            ((model.lb(rxnInd)>=0) && (~any(~ismember(metsL,ActMet))))...% irreversible & reactant side mets all have flux    
           ||...
            ((model.ub(rxnInd)<=0) && (~any(~ismember(metsR,ActMet)))) % irreversible & product side mets all have flux
       
            newLatent = [newLatent;latentCandi(i)];
        end
    end
    % union new latent reaction with previous latent reaction!
    % no new latent rxns?
    if isempty(setdiff(newLatent,latentRxn))
        fprintf('...done! \n');
        break;
    else
        fprintf('...found %d new latent rxns \n',length(setdiff(newLatent,latentRxn)));
        latentRxn = union(latentRxn,newLatent);
    end
    actRxns = union(actRxns,latentRxn);
end
%% maximaze number of "on" latent reaction (like iMAT but without minimize low!)
[A B] = ismember(latentRxn,model.rxns);
latentInd = B(A);
epsilon_f_sorted = epsilon_f(latentInd);
epsilon_r_sorted = epsilon_r(latentInd);
% Creating A matrix
A = [MILProblem.A sparse(size(MILProblem.A,1),2*length(latentRxn));...
    sparse(2*length(latentRxn),size(MILProblem.A,2)) sparse(2*length(latentRxn),2*length(latentRxn))];
for i = 1:length(latentRxn)
    A(i+size(MILProblem.A,1),latentInd(i)) = 1;
    A(i+size(MILProblem.A,1),i+size(MILProblem.A,2)) = model.lb(latentInd(i)) - epsilon_f_sorted(i);
    A(i+size(MILProblem.A,1)+length(latentRxn),latentInd(i)) = 1;
    A(i+size(MILProblem.A,1)+length(latentRxn),i+size(MILProblem.A,2)+length(latentRxn)) = model.ub(latentInd(i)) + epsilon_r_sorted(i);
end
% variable type
vartype1(1:size(MILProblem.A,2),1) = MILProblem.vartype;
vartype2(1:2*length(latentRxn),1) = 'B';
vartype = [vartype1;vartype2];
% Creating csense
csense1(1:size(MILProblem.A,1)) = MILProblem.csense;
csense2(1:length(latentRxn)) = 'G';
csense3(1:length(latentRxn)) = 'L';
csense = [csense1 csense2 csense3];
% Creating lb and ub
lb_y = zeros(2*length(latentRxn),1);
ub_y = ones(2*length(latentRxn),1);
lb = [MILProblem.lb;lb_y];
ub = [MILProblem.ub;ub_y];
% Creating c
c_v = zeros(size(MILProblem.A,2),1);
c_y = ones(2*length(latentRxn),1);
c = [c_v;c_y];
% Creating b
b_s = MILProblem.b;
lb_rh = model.lb(latentInd);
ub_rh = model.ub(latentInd);
b = [b_s;lb_rh;ub_rh];

MILPproblem_latent.A = A;
MILPproblem_latent.b = b;
MILPproblem_latent.c = c;
MILPproblem_latent.lb = lb;
MILPproblem_latent.ub = ub;
MILPproblem_latent.csense = csense;
MILPproblem_latent.vartype = vartype;
MILPproblem_latent.osense = -1;
MILPproblem_latent.x0 = [];
solution = solveCobraMILP_XL(MILPproblem_latent, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 0);
if solution.stat ~= 1
    error('infeasible or violation occured!');
end
fprintf('...total latent fitted: %d \n',sum(solution.int(end-2*length(latentRxn)+1:end)));
Nfit_latent = sum(solution.int(end-2*length(latentRxn)+1:end));
%% minimize total flux
MILPproblem_minFlux = solution2constraint(MILPproblem_latent,solution);
% minimize total flux
% NOTE: this section is specific to the MILP structure in previous
% integration! since we use the V+ variables in the original MILProblem instead of creating new variables
%create a new objective function
% Creating c (objective function)
c_minFlux = zeros(size(MILProblem.A,2),1);
c_minFlux(end-length(Hgenes)-length(model.rxns)+1:end-length(Hgenes)) = ones(length(model.rxns),1);
c = [c_minFlux;zeros(2*length(latentRxn),1)];
MILPproblem_minFlux.c = c;
MILPproblem_minFlux.osense = 1;
solution = solveCobraMILP_XL(MILPproblem_minFlux, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 0);
if solution.stat ~= 1
    error('infeasible or violation occured!');
end
minTotal = solution.obj;
OFD = solution.full(1:length(model.rxns));
end