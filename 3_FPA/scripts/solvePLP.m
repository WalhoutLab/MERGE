function [obj, solution] = solvePLP(model,w, w_label, obj_rxn,default_w, minMax) 
% create a LP for the following problem
% max(myrxn) st w * v <= 1; or the reverse direction
% name it as PLP (Penaltized Linear Problem)
if (nargin < 6)
    minMax = 'max';
end
if (nargin < 5)
    default_w = 0; %the default weight for a reaction if it is not included in any weight input
end

%% step1: match with the weight(w) input
% normally the weight should already be aligned with the rxns in the model,
% but this function also allows input partial weight list with a defined
% label set
% match the weight
w_full = default_w .* ones(length(model.rxns),1);
[A B] = ismember(model.rxns, w_label);
w_full(A) = w(B(A));
%% add the penalty constriants 
model.S(end+1,:) = w_full';
if strcmp(minMax, 'max')
    model.csense(end+1) = 'L';
else
    model.csense(end+1) = 'G';
end
model.b(end+1) = 1; %flux sum is 1 by default; only matters for numeric solubility
model.BiGG(end+1) = {'NA'};
model.metCharges(end+1) = 0;
model.metFormulas(end+1) = {'NA'};
model.metNames(end+1) = {'Nonsense Metabolites which is a constraint'};
model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};
%solve the LP 
model = changeObjective(model,obj_rxn);
solution = optimizeCbModel(model,minMax);
if solution.stat == 1
    obj = solution.obj;
else
    obj = NaN;
end
solution.weight = w_full;
