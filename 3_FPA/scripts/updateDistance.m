function [distMat_new,labels_new] = updateDistance(newDemand,model,distMat,labels,infN)
% this function update the distance matrix to add the distance for the new
% demand rxns (that is for optimizing the producing potential for a
% metabolite)
labels_new = [labels,{[newDemand,'_f']}];
% adding a new dimension to the matrix for the new demand
distMat_new = zeros(1+size(distMat,1),1+size(distMat,2));
distMat_new(1:size(distMat,1),1:size(distMat,2)) = distMat;
% update the distance of rxnX - newDemand 
theMet = model.mets(logical(model.S(:,strcmp(model.rxns,newDemand))));
candidateRxns = setdiff(model.rxns(logical(model.S(strcmp(model.mets,theMet),:))),newDemand);
% add direction label
for i = 1:length(candidateRxns)
    if model.S(strcmp(model.mets,theMet),strcmp(model.rxns,candidateRxns{i})) > 0
        candidateRxns{i} = [candidateRxns{i},'_f'];
    else
        candidateRxns{i} = [candidateRxns{i},'_r'];
    end
end
% remove the irreversible part 
candidateRxns = intersect(candidateRxns,labels);

for i = 1:size(distMat,1)
    distMat_new(i,size(distMat,2)+1) = min(distMat(i,ismember(labels,candidateRxns))) + 1;
    distMat_new(size(distMat,1)+1,i) = infN; %the demand reaction can go to nowhere
end
end