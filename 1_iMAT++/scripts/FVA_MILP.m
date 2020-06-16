function [FVA_lb, FVA_ub] = FVA_MILP(MILPproblem_minFlux, model, targetRxns,parforFlag,BigModel)
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
%    BigModel:          (0 or 1) to indicate if the "big model" mode is
%                       used. This mode is recommanded for all complex models. 
%                       In this mode, we release the MILP strigency 
%                       to gain computational speed. But in general, this mode gives almost
%                       identical flux prediction as the normal mode.

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
if (nargin < 5) 
    BigModel = 0; % by default, do normal FVA
end

if BigModel
    relMipGapTol = 0.001; % we release the MIPgap to 0.1% 
    % this released MipGap only applies to latent step. The strigency of PFD is still kept.
    % users can release the MipGap for PFD manually if needed
else
    relMipGapTol = 1e-12;
end
fprintf('Start to perform the FVA...\n');
% Check if is running on gurobi solver
solverOK = changeCobraSolver('gurobi', 'MILP',0);
if ~solverOK
    fprintf('The solver parameter auto-tuning is not supported for current solver! Please use Gurobi for best performance!\n')
end
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
        if solverOK % for best Gurobi parallel performance, we limit 1 core per parfor thread
            gurobiParameters = struct();
            gurobiParameters.Threads = 1;
            solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
        else % go with default for other solvers
            solution = solveCobraMILP_XL(MILPproblem_minFlux, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
        end
        if solution.stat ~= 1 && solverOK% when failed to solve, we start to tune solver parameter #NOTE: SPECIFIC TO GUROBI SOLVER!%
            gurobiParameters = struct();
            gurobiParameters.Presolve = 0;
            % we dont limit to 1 core to push solver search solution space
            solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
            if solution.stat ~= 1
                gurobiParameters.NumericFocus = 3;
                solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
                if solution.stat ~= 1
                    fprintf('infeasible or violation occured!(%s)\n',targetRxn{:});
                    FVA_lb(i) = nan;
                end
            end
        end        
        if solution.stat == 1
            FVA_lb(i) = solution.obj;
            fprintf('lower boundary of %s found to be %f. \n',targetRxn{:},solution.obj);
        end
        %fprintf('optimizing the the ub of %s...\n',targetRxn{:});
        MILPproblem_minFlux.osense = -1;
        if solverOK % for best Gurobi parallel performance, we limit 1 core per parfor thread
            gurobiParameters = struct();
            gurobiParameters.Threads = 1;
            solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
        else % go with default for other solvers
            solution = solveCobraMILP_XL(MILPproblem_minFlux, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
        end
        if solution.stat ~= 1 && solverOK% when failed to solve, we start to tune solver parameter #NOTE: SPECIFIC TO GUROBI SOLVER!%
            gurobiParameters = struct();
            gurobiParameters.Presolve = 0;
            solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
            if solution.stat ~= 1
                gurobiParameters.NumericFocus = 3;
                solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
                if solution.stat ~= 1
                    fprintf('infeasible or violation occured!(%s)\n',targetRxn{:});
                    FVA_ub(i) = nan;
                end
            end
        end   
        if solution.stat == 1
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
        % when parfor is not used, go with default (use up cores)
        solution = solveCobraMILP_XL(MILPproblem_minFlux, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
        if solution.stat ~= 1 && solverOK% when failed to solve, we start to tune solver parameter #NOTE: SPECIFIC TO GUROBI SOLVER!%
            gurobiParameters = struct();
            gurobiParameters.Presolve = 0;
            solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
            if solution.stat ~= 1
                gurobiParameters.NumericFocus = 3;
                solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
                if solution.stat ~= 1
                    fprintf('infeasible or violation occured!(%s)\n',targetRxn{:});
                    FVA_lb(i) = nan;
                end
            end
        end        
        if solution.stat == 1
            FVA_lb(i) = solution.obj;
            fprintf('lower boundary of %s found to be %f. \n',targetRxn{:},solution.obj);
        end
        %fprintf('optimizing the the ub of %s...\n',targetRxn{:});
        MILPproblem_minFlux.osense = -1;
        solution = solveCobraMILP_XL(MILPproblem_minFlux, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
        if solution.stat ~= 1 && solverOK% when failed to solve, we start to tune solver parameter #NOTE: SPECIFIC TO GUROBI SOLVER!%
            gurobiParameters = struct();
            gurobiParameters.Presolve = 0;
            solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
            if solution.stat ~= 1
                gurobiParameters.NumericFocus = 3;
                solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
                if solution.stat ~= 1
                    fprintf('infeasible or violation occured!(%s)\n',targetRxn{:});
                    FVA_ub(i) = nan;
                end
            end
        end   
        if solution.stat == 1
            FVA_ub(i) = solution.obj;
            fprintf('upper boundary of %s found to be %f. \n',targetRxn{:},solution.obj);
        end
    end
end
end