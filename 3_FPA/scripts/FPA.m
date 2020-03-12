function [fluxEfficiency,fluxEfficiency_plus] = FPA(model,targetRxns,master_expression,distMat,labels,n,manualPenalty,manualDist,maxDist,blockList, constantPenalty,parforFlag,penalty_defined)
%%
if (nargin < 6)
    n = 1.5;
end
if (nargin < 7)
    manualPenalty = {};
end
if (nargin < 8)
    manualDist = {};
end
if (nargin < 9)
    maxDist = max(distMat(~isinf(distMat)));
end
if (nargin < 10)
    blockList = {};
end
if (nargin < 11)
    % constant penaly will not be applied to supercondition 
    constantPenalty = {};
end
if (nargin < 12)
    % allow non-parallel run for metabolite centric optimization;
    % for non-parfor run, the penalty calculation should be pre-defined to
    % save time
    parforFlag = true;
end
if (nargin < 13)
    % allow non-parallel run for metabolite centric optimization;
    % for non-parfor run, the penalty calculation should be pre-defined to
    % save time
    penalty_defined = {};
end
% how to do with no data/partially no data
% how is super condi calculated
%% part1 prepare expression matrix and penalty 
% calculate the penalty from expression data
if parforFlag
    fprintf('Mapping the expression levels to penalties...\n');
    penalty = calculatePenalty(model,master_expression,manualPenalty);
else
    % in non-parfor mode, we need to input penalty (it costs a lot time if
    % repeatedly calculate)
    penalty = penalty_defined;
end

% apply constant penalty to all conditions except for the super condition
% (at the end of the vector) 
[A B] = ismember(model.rxns,constantPenalty(:,1));
penalty(A,1:end-1) = repmat(cell2mat(constantPenalty(B(A),2)),1,size(penalty,2)-1);
%% part2 prepare the distance
% prepare the distance matrix
fprintf('Preparing the distance matrix...\n');
model_irrev = convertToIrreversible(model);
%unify naminclature
tmp_ind = ~cellfun(@(x) any(regexp(x,'_(f|b|r)$')),model_irrev.rxns); %has no suffix
model_irrev.rxns(tmp_ind) = cellfun(@(x) [x,'_f'],model_irrev.rxns(tmp_ind), 'UniformOutput',false);
model_irrev.rxns = regexprep(model_irrev.rxns, '_b$','_r');
%create the full distance matrix
distMat2 = maxDist * ones(length(model_irrev.rxns),length(model_irrev.rxns)); % default distance is the maximum distance
[A B] = ismember(model_irrev.rxns,labels);
%filter Inf 
distMat(distMat > maxDist) = maxDist;
distMat2(A,A) = distMat(B(A),B(A)); %overlay the input distance matrix
if ~isempty(manualDist)
    [A B] = ismember(model_irrev.rxns,manualDist(:,1));
    distMat2(A,:) = repmat(cell2mat(manualDist(B(A),2)),1,size(distMat2,2)); %overlay the manual defined distance (fixed distance reactions) onto the matrix
    distMat2(:,A) = repmat(cell2mat(manualDist(B(A),2)),1,size(distMat2,2))'; %overlay the manual defined distance (fixed distance reactions) onto the matrix
end
% update the distance matrix and label
labels = model_irrev.rxns;
distMat = distMat2;
% reorder expression penalty matrix to match the distance matrix (which
% is reodered to irreversible reactions)
penalty_new = ones(length(labels), size(penalty,2));
for i = 1:length(labels)
    myrxn = labels{i};
    myrxn = myrxn(1:(end-2)); %remove the direction suffix
    penalty_new(i,:) = penalty(strcmp(model.rxns,myrxn),:);
end
%update the penalty matrix to the reordered matrix
penalty = penalty_new;
%% part3 merge penalty and calculate the LP
%form and calculate the LP
fluxEfficiency = cell(length(targetRxns),size(penalty,2));
fluxEfficiency_plus = cell(length(targetRxns),size(penalty,2));
fprintf('FPA calculation started:\n');
if parforFlag
    fprintf(['\n' repmat('.',1,length(targetRxns)) '\n\n']);
    environment = getEnvironment();
    parfor i = 1:length(targetRxns)
        restoreEnvironment(environment);
        myrxn = targetRxns{i};
        doForward = any(strcmp(labels, [myrxn,'_f']));%whether calculate the forward efficiency, according to the distance matrix  
        doReverse = any(strcmp(labels, [myrxn,'_r']));
        efficiencyVector = cell(1,size(penalty,2));
        efficiencyVector_plus = cell(1,size(penalty,2));
        for j = 1:size(penalty,2)
            % block the reactions in the block list
            model_irrev_tmp0 = model_irrev;
            if j < size(penalty,2)
                model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{j})) = 0;
            else %still block when optimizing supertissue // used to not block, now requires a supercond block list (could be empty)
                model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{j})) = 0;
            end
            if doForward 
                model_irrev_tmp = model_irrev_tmp0;
                pDist_f = 1./((1+distMat(strcmp(labels, [myrxn,'_f']),:)).^n);
                w_f = pDist_f' .* penalty(:,j);
                % block the reverse rxns to avoid self-loop
                model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_r'])) = 0;
                [efficiencyVector{1,j}(1),efficiencyVector_plus{1,j}{1}] = solvePLP(model_irrev_tmp,w_f, labels, [myrxn,'_f'],1, 'max');
            else
                efficiencyVector{1,j}(1) = NaN;
                efficiencyVector_plus{1,j}(1) = {''};
            end
            if doReverse 
                model_irrev_tmp = model_irrev_tmp0;
                pDist_r = 1./((1+distMat(strcmp(labels, [myrxn,'_r']),:)).^n);
                w_r = pDist_r' .* penalty(:,j);
                % block the forward rxns to avoid self-loop
                model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_f'])) = 0;
                [efficiencyVector{1,j}(2),efficiencyVector_plus{1,j}{2}] = solvePLP(model_irrev_tmp,w_r, labels, [myrxn,'_r'],1, 'max');
            else
                efficiencyVector{1,j}(2) = NaN;
                efficiencyVector_plus{1,j}(2) = {''};
            end
        end
        fluxEfficiency(i,:) = efficiencyVector;
        fluxEfficiency_plus(i,:) = efficiencyVector_plus;
        fprintf('\b|\n');
    end
else %run the same code on a for loop mode
    for i = 1:length(targetRxns)
        myrxn = targetRxns{i};
        doForward = any(strcmp(labels, [myrxn,'_f']));%whether calculate the forward efficiency, according to the distance matrix  
        doReverse = any(strcmp(labels, [myrxn,'_r']));
        efficiencyVector = cell(1,size(penalty,2));
        efficiencyVector_plus = cell(1,size(penalty,2));
        for j = 1:size(penalty,2)
            % block the reactions in the block list
            model_irrev_tmp0 = model_irrev;
            if j < size(penalty,2)
                model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{j})) = 0;
            else %dont block when optimizing supertissue
            end
            if doForward 
                model_irrev_tmp = model_irrev_tmp0;
                pDist_f = 1./((1+distMat(strcmp(labels, [myrxn,'_f']),:)).^n);
                w_f = pDist_f' .* penalty(:,j);
                % block the reverse rxns to avoid self-loop
                model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_r'])) = 0;
                [efficiencyVector{1,j}(1),efficiencyVector_plus{1,j}{1}] = solvePLP(model_irrev_tmp,w_f, labels, [myrxn,'_f'],1, 'max');
            else
                efficiencyVector{1,j}(1) = NaN;
                efficiencyVector_plus{1,j}(1) = {''};
            end
            if doReverse 
                model_irrev_tmp = model_irrev_tmp0;
                pDist_r = 1./((1+distMat(strcmp(labels, [myrxn,'_r']),:)).^n);
                w_r = pDist_r' .* penalty(:,j);
                % block the forward rxns to avoid self-loop
                model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_f'])) = 0;
                [efficiencyVector{1,j}(2),efficiencyVector_plus{1,j}{2}] = solvePLP(model_irrev_tmp,w_r, labels, [myrxn,'_r'],1, 'max');
            else
                efficiencyVector{1,j}(2) = NaN;
                efficiencyVector_plus{1,j}(2) = {''};
            end
        end
        fluxEfficiency(i,:) = efficiencyVector;
        fluxEfficiency_plus(i,:) = efficiencyVector_plus;
        fprintf('done!\n');
    end
end
