function [solution LPproblem]  = FBA_linearATPm(model,alpha, atpm_name,ATPtargetRxns,osenseStr)
%% set up the objective -- ordinary FBA
if exist('osenseStr', 'var') % Process arguments and set up problem
    if isempty(osenseStr)
        osenseStr = 'max';
    end
else
    if isfield(model, 'osenseStr')
        osenseStr = model.osenseStr;
    else
        osenseStr = 'max';
    end
end
% Figure out objective sense
if strcmpi(osenseStr,'max')
    LPproblem.osense = -1;
elseif strcmpi(osenseStr,'min')
    LPproblem.osense = +1;
else
    error('%s is not a valid osenseStr. Use either ''min'' or ''max''' ,osenseStr);
end
%set up the linear problem
%generate the LHS matrix with new inequality
LPproblem.A = model.S;
%box constraints
LPproblem.lb = model.lb;
LPproblem.ub = model.ub;
%csense
LPproblem.csense = model.csense;
%RHS vector --- assume balance state 
LPproblem.b = zeros(size(LPproblem.A,1),1);
%%
% add the ATPm contraints
LPproblem = addLinearATPm(LPproblem, alpha, model,atpm_name,ATPtargetRxns);
%% linear objective coefficient 
c_atp = zeros(length(model.rxns),1);
LPproblem.c = [model.c;c_atp];
% complie csense
LPproblem.csense = columnVector(LPproblem.csense);
%Double check that all inputs are valid:
if ~(verifyCobraProblem(LPproblem, [], [], false) == 1)
    warning('invalid problem');
    return;
end

solution = solveCobraLP_XL(LPproblem);
end
