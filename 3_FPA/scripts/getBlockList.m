function blockList = getBlockList(model,FVAtbl)
%% To be completed!!

model_irrev = convertToIrreversible(model);
%unify naminclature
tmp_ind = ~cellfun(@(x) any(regexp(x,'_(f|b|r)$')),model_irrev.rxns); %has no suffix
model_irrev.rxns(tmp_ind) = cellfun(@(x) [x,'_f'],model_irrev.rxns(tmp_ind), 'UniformOutput',false);
model_irrev.rxns = regexprep(model_irrev.rxns, '_b$','_r');
X_rxns = model_irrev.rxns(contains(model_irrev.rxns,'_X_'));
for i = 1:length(tissueLabel)
    blockList{i} = confidenceTbl.Properties.RowNames((confidenceTbl.(tissueLabel{i}) <= -2));
    if strcmp(tissueLabel{i},'Intestine')
        blockList{i} = [blockList{i};X_rxns]; %block all X tissue reactions for intestine
    end     
end
blockList{end+1} = {};%for supercond, dont block anything when optimizing X tissue
blockList_I = blockList;
blockList_I{end} = X_rxns;%s
