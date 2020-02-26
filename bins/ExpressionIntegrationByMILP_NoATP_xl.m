function [OpenGene, OpenedHReaction,ClosedLReaction,solution, MILPproblem] = ExpressionIntegrationByMILP_NoATP_xl(model, RHNames, RLNames,HGeneList,epsilon_f, epsilon_r, tol, core, logfile, runtime)
% the script is based on iMAT function with significant alteration.
%the overall strategy in this algrithom is to formulate the gene expression
%data into a MILP where number of high expressed reaction(genes) is
%maximized and low expressed reaction(genes) is minimized.
%
% USAGE:
%
%    tissueModel = iMAT(model, expressionRxns, threshold_lb, threshold_ub)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    RHindex:           cell with reaction names that are labled as high expression; If in gene maximization algorithm, the reactions in this list will be associated with genes.
%    RLNames:           cell with reaction names that are labled as low
%                       expression; actually is zero level in the
%                       quadra-level catagory and needs to be minimized
%
%
% OPTIONAL INPUTS:
%    tol:               minimum flux threshold for "expressed" reactions
%                       (default 1e-8)
%    core:              cell with reaction names (strings) that are manually put in
%                       the high confidence set (default - no core reactions)
%    logfile:           name of the file to save the MILP log (string)
%    runtime:           maximum solve time for the MILP (default value - 7200s)
%    epsilon_f/r:           a vector for all the high flux reactions, need to
%                       be in the same length of RHnames and in the same
%                       order; epsilons are in absolute value
% OUTPUT:
%    tissueModel:       extracted model
%
% `Zur et al. (2010). iMAT: an integrative metabolic analysis tool. Bioinformatics 26, 3140-3142.`
%
% .. Author: - Implementation adapted from the cobra toolbox
% (createTissueSpecificModel.m) by S. Opdam and A. Richelle, May 2017
% Xuhang Li, Dec 2018 

if nargin < 13 || isempty(runtime)
    runtime = 7200;
    %runtime = 60;
end
if nargin < 12 || isempty(logfile)
    logfile = 'MILPlog';
end

if nargin < 11 || isempty(core)
    core={};
end

if nargin < 10 || isempty(tol)
    tol = 1e-8;
end
%% step1: formulate the high expression index(RHindex) and low expression
%index(RLindex) 
RHindex = find(ismember(model.rxns, RHNames));
RLindex = find(ismember(model.rxns, RLNames));
epsilon_f_sorted = epsilon_f(ismember(model.rxns,RHNames));
epsilon_r_sorted = epsilon_r(ismember(model.rxns,RHNames));


%Manually add defined core reactions to the core
if ~isempty(core)
    for i = 1:length(core)
        rloc = find(ismember(model.rxns, core{i}));
        if ~isempty(rloc) && isempty(intersect(RHindex,rloc))
            RHindex(end+1) = rloc;
        end
        %core set has higher priority than low expression label
        %to overwrite possible core reactions that are labled as low
        %expressed 
        if ~isempty(rloc) && ~isempty(intersect(RLindex,rloc))
            RLindex(RLindex == intersect(RLindex,rloc)) = [];
        end
        if isempty(rloc)
            disp(['Manual added core reaction: ', core{i}, ' not found'])
        end
    end
end
%now the expression and core list inputs are integrated into RHindex
%and RLindex
%% step2: formulate the MILP
%start to generate the MILP input
S = model.S;
lb = model.lb;
ub = model.ub;

% formulate the basic MILP with objective function undefined 
% -------------------------
% Creating A matrix
A = sparse(size(S,1)+2*length(RHindex)+2*length(RLindex),size(S,2)+2*length(RHindex)+length(RLindex));
[m,n,s] = find(S);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end

for i = 1:length(RHindex)
    A(i+size(S,1),RHindex(i)) = 1;
    A(i+size(S,1),i+size(S,2)) = lb(RHindex(i)) - epsilon_f_sorted(i);
    A(i+size(S,1)+length(RHindex),RHindex(i)) = 1;
    A(i+size(S,1)+length(RHindex),i+size(S,2)+length(RHindex)+length(RLindex)) = ub(RHindex(i)) + epsilon_r_sorted(i);
end

for i = 1:length(RLindex)
    A(i+size(S,1)+2*length(RHindex),RLindex(i)) = 1;
    A(i+size(S,1)+2*length(RHindex),i+size(S,2)+length(RHindex)) = lb(RLindex(i));
    A(i+size(S,1)+2*length(RHindex)+length(RLindex),RLindex(i)) = 1;
    A(i+size(S,1)+2*length(RHindex)+length(RLindex),i+size(S,2)+length(RHindex)) = ub(RLindex(i));
end
% Creating csense
csense1(1:size(S,1)) = model.csense;
csense2(1:length(RHindex)) = 'G';
csense3(1:length(RHindex)) = 'L';
csense4(1:length(RLindex)) = 'G';
csense5(1:length(RLindex)) = 'L';
csense = [csense1 csense2 csense3 csense4 csense5];
% Creating lb and ub
lb_y = zeros(2*length(RHindex)+length(RLindex),1);
ub_y = ones(2*length(RHindex)+length(RLindex),1);
lb = [lb;lb_y];
ub = [ub;ub_y];
% Creating b
b_s = model.b;
lb_rh = lb(RHindex);
ub_rh = ub(RHindex);
lb_rl = lb(RLindex);
ub_rl = ub(RLindex);
b = [b_s;lb_rh;ub_rh;lb_rl;ub_rl];
% Creating vartype
vartype1(1:size(S,2),1) = 'C';
vartype2(1:2*length(RHindex)+length(RLindex),1) = 'B';
vartype = [vartype1;vartype2];
% create MILP
MILPproblem.A = A;
MILPproblem.b = b;
MILPproblem.lb = lb;
MILPproblem.ub = ub;
MILPproblem.csense = csense;
MILPproblem.vartype = vartype;
MILPproblem.osense = -1;
MILPproblem.x0 = [];
% ---------------------------
%%
% add the ATPm contraints
%MILPproblem = addLinearATPm(MILPproblem, alpha, model,atpm_name,ATPtargetRxns);

% add the gene expression maximization constraint
MILPproblem = addGeneMaximization(MILPproblem, HGeneList,RHNames,RLNames, model);

% Creating c (objective function)
c_v = zeros(size(S,2),1);
c_yh1 = zeros(length(RHindex),1);
c_yl = ones(length(RLindex),1);
c_yh2 = zeros(length(RHindex),1);
c_gene = ones(length(HGeneList),1);
c = [c_v;c_yh1;c_yl;c_yh2;c_gene];
MILPproblem.c = c;
%% Step 3 generate the output
solution = solveCobraMILP_XL(MILPproblem, 'timeLimit', runtime, 'logFile', logfile, 'printLevel', 0);
if solution.stat ~= 1
    error('infeasible or violation occured!');
end
fprintf('...total high gene fitted: %d \n',sum(solution.int(end-length(c_gene)+1:end)));
yH = boolean(solution.int(1:length(RHindex)) + solution.int((length(RHindex)+length(RLindex)+1):(2*length(RHindex)+length(RLindex))));
nameH = model.rxns(ismember(model.rxns,RHNames));
OpenedHReaction = nameH(yH);
nameL = model.rxns(ismember(model.rxns,RLNames));
yL = boolean(solution.int((length(RHindex)+1):(length(RHindex)+length(RLindex))));
ClosedLReaction = nameL(yL);
OpenGene = HGeneList(boolean(solution.int((2*length(RHindex)+length(RLindex))+1:end)));
%x = solution.cont;
%rxnRemList = model.rxns(abs(x) < tol);
end

