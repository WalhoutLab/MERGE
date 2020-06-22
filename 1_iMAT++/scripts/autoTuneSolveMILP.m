function [solution] = autoTuneSolveMILP(MILP,solverOK,relMipGapTol,marker, limitCore)
% The automatic solver optimizer for parallel MILP calculation in FAA

% double the default time limit (we try to get FVA done for as many reaction as
% possible)

% non-optimal solution returns NaN

if (nargin < 5) 
    limitCore = 1; % by default, we assume parfor is used and we limit the solver core usage
end

if solverOK && limitCore% for best Gurobi parallel performance, we limit 1 core per parfor thread
    gurobiParameters = struct();
    gurobiParameters.Threads = 1;
    solution = solveCobraMILP_XL(MILP, gurobiParameters,'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
else % go with default for other solvers
    solution = solveCobraMILP_XL(MILP, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
end

if solution.stat ==0 && solverOK% when failed to solve (infeasible), we start to tune solver parameter #NOTE: SPECIFIC TO GUROBI SOLVER!%
    gurobiParameters = struct();
    gurobiParameters.Presolve = 0;
    % we dont limit to 1 core to push solver search solution space
    solution = solveCobraMILP_XL(MILP,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
    if solution.stat == 0 % still infeasible
        gurobiParameters.NumericFocus = 3;
        solution = solveCobraMILP_XL(MILP,gurobiParameters, 'timeLimit', 1200, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
    end
end

if solution.stat ~= 1 % give up and report error
    fprintf('infeasible or unbounded occured!(%s) (error code: %d)\n',marker,solution.stat);
    solution.obj = nan;
end
end