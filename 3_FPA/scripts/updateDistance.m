function [distMat_new,labels_new] = updateDistance(newDemand,model,distMat,labels,infN)
% update the distance matrix for the newly added demand reaction. Used in
% the metabolite-centric FPA
%
% USAGE:
%
%    [distMat_new,labels_new] = updateDistance(newDemand,model,distMat,labels,infN)
%
% INPUT:
%    newDemand:           newly added demand reaction (reaction name) 
%    model:               the COBRA model
%    distMat:             original distance matrix [note: need to be the
%                         raw distance matrix (from rxn A to rxn B,
%                         direction matters)]
%    labels:              labels (the reaction name) of the original
%                         distance matrix
%    infN:                maximum distance (representing infinity)
% OUTPUT:
%    distMat_new:         the new distance matrix 
%    labels_new:          labels for the new distance matrix
%
% .. Author: Xuhang Li, Mar 2020 

% add the new demand reaction to the end of the distance matrix and label
labels_new = [labels,{[newDemand,'_f']}];
% adding a new dimension to the matrix for the new demand
distMat_new = zeros(1+size(distMat,1),1+size(distMat,2));
distMat_new(1:size(distMat,1),1:size(distMat,2)) = distMat;
% update the distance of the newDemand [which is min(drained met to a
% reaction) + 1]
% first find all reaction that contains the drained metabolite
theMet = model.mets(logical(model.S(:,strcmp(model.rxns,newDemand))));
candidateRxns = setdiff(model.rxns(logical(model.S(strcmp(model.mets,theMet),:))),newDemand);
% add direction label
for i = 1:length(candidateRxns)
    if model.S(strcmp(model.mets,theMet),strcmp(model.rxns,candidateRxns{i})) > 0% only the producing direction is considered
        candidateRxns{i} = [candidateRxns{i},'_f'];
    else
        candidateRxns{i} = [candidateRxns{i},'_r'];
    end
end
% remove the irreversible part (some reactions are actully irreversible)
candidateRxns = intersect(candidateRxns,labels);
% take the minimum and make output
for i = 1:size(distMat,1)
    distMat_new(i,size(distMat,2)+1) = min(distMat(i,ismember(labels,candidateRxns))) + 1;
    distMat_new(size(distMat,1)+1,i) = infN; % the demand reaction can go to nowhere
end
end