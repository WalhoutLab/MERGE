function [newModel] = collapseX(model, compID, fluxDistribution)
% this function is to collapse one compartment in a conjoined model to get
% a collapsed model; the collapsed compartment will be represented as set
% of psudo-exchange reations that are in place of the original transport
% reaction; these exchange reaction phenocopy the original demand of each
% metabolite in the collapsed compartment 
newModel = model;
%% collapse metabolite and reactions
Xrxns = model.rxns(cellfun(@(x) ~isempty(regexp(x, ['_',compID,'$'], 'once')),model.rxns));
newModel = removeRxns(newModel,Xrxns);
Xmets = model.mets(cellfun(@(x) ~isempty(regexp(x, ['_',compID,'\[.\]$'], 'once')),model.mets));
newModel = removeMetabolites(newModel,Xmets);
%collapsing genes in X is skipped since reshaping the GR matrix is slow
%Xgenes = cellfun(@(x) ~isempty(regexp(x, ['_',compID,'$'], 'once')),model.genes);
%newModel =  removeFieldEntriesForType(newModel,Xgenes,'genes',numel(newModel.genes));


%note: a scar is left: the sideConstraints (the NonsenseMets) for X
%compartement is not deleted. but they are now essentially all zero vector,
%which will not influence the modeling 
%% add psudo-exchange rxns
% calculate the flux demand for each metabolite based on transport
% reactions
% find the X_E metabolites
numericalTol = 0; % tolerance could be larger or equal to solver tolerance; zero gives strigency but may cause solver instability
allExMets = model.mets(cellfun(@(x) ~isempty(regexp(x, '\[e\]$', 'once')),model.mets));
allXmets = model.mets(any(model.S(:,ismember(model.rxns,Xrxns)),2)); % all metablites used in a X reaction
targetMets = intersect(allExMets, allXmets);
% add reactions
for i = 1:length(targetMets)
    % determine the lb and ub
    % consider the total comsuming flux
    associatedRxns = model.rxns(any(model.S(strcmp(model.mets, targetMets{i}),:),1));
    associatedRxns = intersect(associatedRxns,Xrxns);
    demandFlux = -(model.S(strcmp(model.mets, targetMets{i}),ismember(model.rxns,associatedRxns)) * fluxDistribution(ismember(model.rxns,associatedRxns)));
    if demandFlux == 0
        deltaV = 0;
    else
        deltaV = numericalTol;
    end
    % construct the reaction
    newModel = addReaction(newModel,['EX',model.MetMachineID{strcmp(targetMets{i},model.mets)}(2:end),['_',compID]],...
        'metaboliteList',targetMets(i),...
        'stoichCoeffList',-1,...
        'reversible',1,...
        'subSystem', 'Collapsed compartment X',...
        'lowerBound',demandFlux-deltaV,...
        'upperBound',demandFlux+deltaV,'printLevel',0);
end  

end