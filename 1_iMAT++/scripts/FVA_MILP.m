function [FVA_lb, FVA_ub] = FVA_MILP(MILPproblem_minFlux, model, targetRxns,parforFlag)
% perform FVA given a setup MILProblem and target reactions. The first nRxn
% variables in the MILProblem must be the same as reactions in the input
% model.
%
% USAGE:
%
%    [lb, ub] = FVA_MILP(MILProblem, model, targetRxns,parforFlag)
%
% INPUTS:
%    MILPproblem_minFlux: the input MILP problem (COBRA MILP structure).The
%                       MILP should be readily constrained for FVA
%                       calculation. For example, the total flux cap should
%                       be already set.
%    model:             input model (COBRA model structure)
%    targetRxns:        cell of target reactions to perform FVA on
%    parforFlag:        (0 or 1) whether to use parallel computing
%
% OUTPUT:
%   ub:                 a vector of upper boundaries of queried reactions
%   lb:                 a vector of lower boundaries of queried reactions
%
% Additional Notice:    Please make sure the S matrix of the input MILP follows the structure of iMAT++ MILP. Some variables such as absolute flux proxy will be assumed to be at specifc positions, so errors will occur if the S matrix is not formed as standard iMAT++. 
%
% `Yilmaz et al. (2020). Final Tittle and journal.
% .. Author: - Xuhang Li, Mar 2020
if nargin < 3 || isempty(targetRxns)
    targetRxns = model.rxns;
end
if nargin < 4 || isempty(parforFlag)
    parforFlag = true;
end
fprintf('Start to perform the FVA...\n');
%% analyze FVA
if parforFlag
    environment = getEnvironment();
    MILPproblem_minFlux_ori = MILPproblem_minFlux;
    parfor i = 1:length(targetRxns)
        restoreEnvironment(environment);
        MILPproblem_minFlux = MILPproblem_minFlux_ori;
        targetRxn = targetRxns(i);
        FluxObj = find(ismember(model.rxns,targetRxn)); 
        %create a new objective function
        c = zeros(size(MILPproblem_minFlux.A,2),1);
        c(FluxObj) = 1;
        MILPproblem_minFlux.c = c;
        MILPproblem_minFlux.osense = 1;
        %fprintf('optimizing for the lb of %s...\n',targetRxn{:});
        solution = solveCobraMILP_XL(MILPproblem_minFlux, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 0);
        if solution.stat ~= 1
            error('infeasible or violation occured!');
        else
            FVA_lb(i) = solution.obj;
            fprintf('lower boundary of %s found to be %f. \n',targetRxn{:},solution.obj);
        end
        %fprintf('optimizing the the ub of %s...\n',targetRxn{:});
        MILPproblem_minFlux.osense = -1;
            solution = solveCobraMILP_XL(MILPproblem_minFlux, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 0);
        if solution.stat ~= 1
            error('infeasible or violation occured!');
        else
            FVA_ub(i) = solution.obj;
            fprintf('upper boundary of %s found to be %f. \n',targetRxn{:},solution.obj);
        end
    end
else %same thing but in for loop
    for i = 1:length(targetRxns)
        targetRxn = targetRxns(i);
        FluxObj = find(ismember(model.rxns,targetRxn)); 
        %create a new objective function
        c = zeros(size(MILPproblem_minFlux.A,2),1);
        c(FluxObj) = 1;
        MILPproblem_minFlux.c = c;
        MILPproblem_minFlux.osense = 1;
        %fprintf('optimizing for the lb of %s...\n',targetRxn{:});
        solution = solveCobraMILP_XL(MILPproblem_minFlux, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 0);
        if solution.stat ~= 1
            error('infeasible or violation occured!');
        else
            FVA_lb(i) = solution.obj;
            fprintf('lower boundary of %s found to be %f. \n',targetRxn{:},solution.obj);
        end
        %fprintf('optimizing the the ub of %s...\n',targetRxn{:});
        MILPproblem_minFlux.osense = -1;
            solution = solveCobraMILP_XL(MILPproblem_minFlux, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 0);
        if solution.stat ~= 1
            error('infeasible or violation occured!');
        else
            FVA_ub(i) = solution.obj;
            fprintf('upper boundary of %s found to be %f. \n',targetRxn{:}, solution.obj);
        end
    end
end
end