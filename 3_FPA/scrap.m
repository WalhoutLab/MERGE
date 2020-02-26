targetRxns = [targetRxns; alltarget(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE|SNK|DMN)','once')),alltarget));targetMets];
%%
compare = table();
compare.safak = nan(size(FPAtbl,1),1);
compare.hang = nan(size(FPAtbl,1),1);
for j = 1:length(tissueLabel)
    tissue = tissueLabel{j};
    for i = 1:size(FPAtbl,1)
        compare.safak(i) = FPAtbl.(tissue)(i);
        if strcmp(FPAtbl.Properties.RowNames{i}(end),'f')
            compare.hang(i) = relFP_f(strcmp(targetRxns,FPAtbl.Properties.RowNames{i}(1:end-1)),strcmp(tissueLabel,tissue));
        else
            compare.hang(i) = relFP_r(strcmp(targetRxns,FPAtbl.Properties.RowNames{i}(1:end-1)),strcmp(tissueLabel,tissue));
        end
    end
    compare.hang(isnan(compare.hang)) = 0;
    compare.safak(isnan(compare.safak)) = 0;
    fprintf('for tissue %s find %d/%d identical values;\n',tissue,sum(abs(compare.hang - compare.safak)<=0.0001),size(FPAtbl,1));
    wrongrxn{j} = FPAtbl.Properties.RowNames(~(abs(compare.hang - compare.safak)<=0.0001));
end
%%
check = intersect(intersect(intersect(intersect(intersect(intersect(wrongrxn{1},wrongrxn{2}),wrongrxn{3}),wrongrxn{4}),wrongrxn{5}),wrongrxn{6}),wrongrxn{7})