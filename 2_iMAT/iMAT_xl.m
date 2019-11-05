function [solution, MILPproblem,Hreactions,Lreactions] = iMAT_xl(model, expressionRxns,epsilon_f,epsilon_r, threshold_lb, threshold_ub, core, logfile, runtime)
% this algorithm is adapted from original iMAT as following: 
%--------------
% Uses the iMAT algorithm (`Zur et al., 2010`) to extract a context
% specific model using data. iMAT algorithm find the optimal trade-off
% between inluding high-expression reactions and removing low-expression reactions.
%
% USAGE:
%
%    tissueModel = iMAT(model, expressionRxns, threshold_lb, threshold_ub)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    expressionRxns:    reaction expression, expression data corresponding to model.rxns.
%                       Note : If no gene-expression data are
%                       available for the reactions, set the
%                       value to -1
%    threshold_lb:      lower bound of expression threshold, reactions with
%                       expression below this value are "non-expressed"
%    threshold_ub:      upper bound of expression threshold, reactions with
%                       expression above this value are "expressed"
%
%
% OPTIONAL INPUTS:
%    tol:               minimum flux threshold for "expressed" reactions
%                       (default 1e-8)
%    core:              cell with reaction names (strings) that are manually put in
%                       the high confidence set (default - no core reactions)
%    logfile:           name of the file to save the MILP log (string)
%    runtime:           maximum solve time for the MILP (default value - 7200s)
%
% OUTPUT:
%    tissueModel:       extracted model
%
% `Zur et al. (2010). iMAT: an integrative metabolic analysis tool. Bioinformatics 26, 3140-3142.`
%
% .. Author: - Implementation adapted from the cobra toolbox
% (createTissueSpecificModel.m) by S. Opdam and A. Richelle, May 2017
%--------------
% the modification enables dynamic epsilon for different rxns. modified
% lines are commented
% Xuhang Li, 07/30/2019

if nargin < 9 || isempty(runtime)
    runtime = 7200;
    %runtime = 60;
end
if nargin < 8 || isempty(logfile)
    logfile = 'MILPlog';
end
if nargin < 7 || isempty(core)
    core={};
end

%sort epsilon vector according to reaction order

RHindex = find(expressionRxns >= threshold_ub);
Hreactions = model.rxns(RHindex);
RLindex = find(expressionRxns >= 0 & expressionRxns < threshold_lb);
Lreactions = model.rxns(RLindex);
epsilon_f_sorted = epsilon_f(RHindex);
epsilon_r_sorted = epsilon_r(RHindex);
%Manually add defined core reactions to the core
%core reaction function is preserved even it is not used here

if ~isempty(core)
    for i = 1:length(core)
        rloc = find(ismember(model.rxns, core{i}));
        if ~isempty(rloc) && isempty(intersect(RHindex,rloc))
            RHindex(end+1) = rloc;
        end
        %to overwrite unrealistic RLindex caused by artificial gene exp
        %settings!!!!! specific to ROS project!!!!!!
        if ~isempty(rloc) && ~isempty(intersect(RLindex,rloc))
            RLindex(RLindex == intersect(RLindex,rloc)) = [];
        end
        if isempty(rloc)
            disp(['Manual added core reaction: ', core{i}, ' not found'])
        end
    end
end

S = model.S;
lb = model.lb;
ub = model.ub;

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
    A(i+size(S,1)+length(RHindex),i+size(S,2)+length(RHindex)+length(RLindex)) = ub(RHindex(i)) + epsilon_r_sorted(i);%dynamic epsilon is used
end

for i = 1:length(RLindex)
    A(i+size(S,1)+2*length(RHindex),RLindex(i)) = 1;
    A(i+size(S,1)+2*length(RHindex),i+size(S,2)+length(RHindex)) = lb(RLindex(i));
    A(i+size(S,1)+2*length(RHindex)+length(RLindex),RLindex(i)) = 1;
    A(i+size(S,1)+2*length(RHindex)+length(RLindex),i+size(S,2)+length(RHindex)) = ub(RLindex(i));
end

% Creating csense
csense1(1:size(S,1),1) = model.csense;
csense2(1:length(RHindex),1) = 'G';
csense3(1:length(RHindex),1) = 'L';
csense4(1:length(RLindex),1) = 'G';
csense5(1:length(RLindex),1) = 'L';
csense = [csense1;csense2;csense3;csense4;csense5];

% Creating lb and ub
lb_y = zeros(2*length(RHindex)+length(RLindex),1);
ub_y = ones(2*length(RHindex)+length(RLindex),1);
lb = [lb;lb_y];
ub = [ub;ub_y];

% Creating c
c_v = zeros(size(S,2),1);
c_y = ones(2*length(RHindex)+length(RLindex),1);
c = [c_v;c_y];

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

MILPproblem.A = A;
MILPproblem.b = b;
MILPproblem.c = c;
MILPproblem.lb = lb;
MILPproblem.ub = ub;
MILPproblem.csense = csense;
MILPproblem.vartype = vartype;
MILPproblem.osense = -1;
MILPproblem.x0 = [];

solution = solveCobraMILP_XL(MILPproblem, 'timeLimit', runtime, 'logFile', logfile, 'printLevel', 0);
    
end

