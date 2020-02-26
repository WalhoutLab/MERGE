function penalty = calculatePenalty(model,master_expression,manualPenalty)
% impose expression to GPR with AND blocks
% in the current mapping algorithm, the no data is labeled as NaN, and
% partial data is set as the no data gene ignored 
levels = cell(length(model.rxns),length(master_expression));
status = -1*ones(length(model.rxns),length(master_expression));
for i = 1:length(master_expression)
    expression = master_expression{i};
    expression.value = expression.value + 1; %add 1 pseudocount
    [levels(:,i), status(:,i)] = gene_to_reaction_levels(model, expression.genes, expression.value, @min, @(x,y)(x+y));
end
% calculate penalty
% the NaNs are set as default penalty (1)
normalizedLevel = nan(size(levels,1),size(levels,2));
penalty = ones(size(levels,1),size(levels,2));%initial with default penalty
for i = 1:length(model.rxns)
    lenL = zeros(length(master_expression));
    for j=1:length(master_expression)
        lenL(j)=length(levels{i,j});
    end
    if any(lenL ~= lenL(1)) %some GPR parsing error
        error('GPR length unequal');
    else
        %normalize every block
        stackM = nan(lenL(1),length(master_expression)); %stack in a matrix
        for z = 1:lenL(1)
            for s = 1:length(master_expression)
                stackM(z,s) = levels{i,s}(z);
            end
        end
        %normalize the matrix
        stackM = stackM ./ max(stackM,[],2);
        %the safak-supercond equals to max(stackM,[],2). so it will always
        %be 1
        %pick the minimal value as expression for each sample
        normalizedLevel(i,:) = min(stackM,[],1);
        if ~any(isnan(normalizedLevel(i,:)))
            penalty(i,:) = 1./normalizedLevel(i,:);
        elseif any(isnan(normalizedLevel(i,:))) && ~all(isnan(normalizedLevel(i,:)))
            error('partial None Penalty for some conditions, check!');
        end
    end
end

%add a super condition column in the penalty matrix
supercond = ones(size(penalty,1),1);
penalty = [penalty,supercond];

%some reactions like exchange reactions are non-GPR reaction; this could be
%manually set as needed
if exist('manualPenalty','var')
    for i = 1:size(manualPenalty,1)
        penalty(strcmp(model.rxns,manualPenalty{i,1}),:) = manualPenalty{i,2};
    end
end
end