function solution  = FBA_dynamic_biomass(model,biomass_drain,biomass_reference, lRange,uRange,osenseStr)
%% set up the objective
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
biomassPartInd = [find(cellfun(@(x) ~isempty(regexp(x, 'biomass_L_part_')),model.rxns));find(cellfun(@(x) ~isempty(regexp(x, 'biomass_R_part_')),model.rxns))];
biomassDrainInd = find(strcmp(model.rxns,biomass_drain));
biomassPartName = {};
for i = biomassPartInd'
    str1 = strsplit(model.rxns{i},'part_');
    biomassPartName = [biomassPartName,str1(2)];
end
MetInd = cellfun(@(x) find(strcmp(model.mets, x)),biomassPartName);
coef = abs(model.S(MetInd,strcmp(model.rxns,biomass_reference)));
LPproblem.A = [LPproblem.A;sparse(2*length(coef),size(LPproblem.A,2))];
for i = 1:length(coef)
    LPproblem.A(size(model.S,1)+i,biomassPartInd(i)) = -1;
    LPproblem.A(size(model.S,1)+i,biomassDrainInd) = coef(i) * uRange;
end
for i = 1:length(coef)
    LPproblem.A(size(model.S,1)+length(coef)+i,biomassPartInd(i)) = -1;
    LPproblem.A(size(model.S,1)+length(coef)+i,biomassDrainInd) = coef(i) * lRange;
end
%linear objective coefficient
LPproblem.c = model.c;
%box constraints
LPproblem.lb = model.lb;
LPproblem.ub = model.ub;
%csense
LPproblem.csense=[model.csense,repmat('G',1,length(coef)),repmat('L',1,length(coef))];
LPproblem.csense = columnVector(LPproblem.csense);
%RHS vector --- assume balance state 
LPproblem.b = zeros(size(LPproblem.A,1),1);

%Double check that all inputs are valid:
if ~(verifyCobraProblem(LPproblem, [], [], false) == 1)
    warning('invalid problem');
    return;
end

solution = solveCobraLP_XL(LPproblem);
end
