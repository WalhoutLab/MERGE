function [solutionDel, solStatus] = linearMOMA_xl(modelDel, V_wt)
% Performs a linear version of the MOMA (minimization of metabolic adjustment) approach
%
% USAGE:
%
%    [solutionDel, solutionWT, totalFluxDiff, solStatus] = linearMOMA(modelWT, modelDel, osenseStr, minFluxFlag, verbFlab)
%
% INPUTS:
%    modelWT:          Wild type model
%    modelDel:         Deletion strain model
%
% OPTIONAL INPUTS:
%    osenseStr:        Maximize ('max') / minimize ('min') (Default = 'max')
%    minFluxFlag:      Minimize the absolute value of fluxes in the optimal MOMA
%                      solution (Default = false)
%    verbFlag:         Verbose output (Default = false)
%
% OUTPUTS:
%    solutionDel:      Deletion solution structure
%    solutionWT:       Wild-type solution structure
%    totalFluxDiff:    Value of the linear MOMA objective, i.e. :math:`\sum |v_{wt}-v_{del}|`
%    solStatus:        Solution status - solves the problem: (`f_wt` is the optimal wild type objective value found by FBA)
%
% .. math::
%     min ~&~  \sum |v_{wt} - v_{del}| \\
%         ~&~ S_{wt}v_{wt} = 0 \\
%         ~&~ lb_{wt} \leq v_{wt} \leq ub_{wt} \\
%         ~&~ c_{wt}^T v_{wt} = f_{wt} \\
%         ~&~ S_{del}v_{del} = 0 \\
%         ~&~ lb_{del} \leq v_{del} \leq ub_{del}
%
% NOTE:
%
%    1) This formulation allows for selecting the most appropriate
%    optimal wild type FBA solution as the starting point as opposed to
%    picking an arbitrary starting point (original MOMA implementation).
%
%    2) The reaction sets in the two models do not have to be equal as long as
%    there is at least one reaction in common
%
% .. Author: - Markus Herrgard 11/7/06

% LP solution tolerance

[nMets2,nRxns2] = size(modelDel.S);
solutionDel.f = [];
solutionDel.x = [];
solutionDel.stat = -1;

% Variables in the following problem are
% x = [v1;v2;delta+;delta-]
% where v1 = wild type flux vector
%       v2 = deletion strain flux vector
%       delta+ = v1 - v2
%       delta- = v2 - v1


% Construct the LHS matrix
% Rows:
% 1: Swt*v1 = 0 for the wild type
% 2: Sdel*v2 = 0 for the deletion strain
% 3: delta+ >= v1-v2
% 4: delta- >= v2-v1
% 5: c'v1 = f1 (wild type)
A =  [modelDel.S sparse(nMets2,2*nRxns2);
     eye(nRxns2) eye(nRxns2) sparse(nRxns2,nRxns2);
     -eye(nRxns2) sparse(nRxns2,nRxns2) eye(nRxns2)];
% Construct the RHS vector
b = [zeros(nMets2,1);V_wt; -V_wt];
% Construct the objective (sum of all delta+ and delta-)
c = [zeros(nRxns2,1);ones(2*nRxns2,1)];
% Construct the ub/lb
% delta+ and delta- are in [0 10000]
lb = [modelDel.lb;zeros(2*nRxns2,1)];
ub = [modelDel.ub;10000*ones(2*nRxns2,1)];
% Construct the constraint direction vector (G for delta's, E for
% everything else)
csense(1:(nMets2)) = 'E';
csense((nMets2)+1:(nMets2+2*nRxns2)) = 'G';
% Solve the linearMOMA problem
[LPproblem.A,LPproblem.b,LPproblem.c,LPproblem.lb,LPproblem.ub,LPproblem.csense,LPproblem.osense] = deal(A,b,c,lb,ub,csense,1);
LPsolution = solveCobraLP(LPproblem);

if (LPsolution.stat > 0)
    solutionDel.x = LPsolution.full(1:nRxns2);
    solutionDel.stat = LPsolution.stat;
    solutionDel.f = LPsolution.obj;
end
solStatus = LPsolution.stat;
