distance = readtable('distance_rxn2rxn.tsv','FileType','text');
%% convert the file to a matrix
bak = distance;
distance = table2cell(distance);
labels = unique(distance(:,1));
distance_num = cell2mat(distance(:,3));
distMat = nan(length(labels),length(labels));
for i = 1:length(labels)
    tmp_ind = strcmp(distance(:,1),labels{i});
    tmp_label = distance(tmp_ind,2);
    tmp_num = distance_num(tmp_ind);
    [A B] = ismember(labels, tmp_label);
    distMat(i,A) = tmp_num(B(A));
end
for i = 1:length(labels)-1
    for j = i+1:length(labels)
        if isnan(distMat(i,j))
            error('!')
        end
    end
end
%reformat label
labels = cellfun(@(x) [x(1:end-1),'_',x(end)],labels,'UniformOutput',false);
% expand the matrix
for i = 1:length(labels)
    for j = 1:i-1
        distMat(i,j) = distMat(j,i);
    end
    distMat(i,i) = 0;
end
save('distanceMatrix.mat','distMat','labels');
