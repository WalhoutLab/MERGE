function mytbl = listRxn(model,flux,myMet)
myrxns = model.rxns(any(model.S(strcmp(model.mets,myMet),:),1));
myInds = any(model.S(strcmp(model.mets,myMet),:),1);
myfluxP = flux(ismember(model.rxns,myrxns)) .*  model.S(strcmp(model.mets,myMet),myInds)';
for i = 1:length(myrxns)
    mytbl(i,1) = myrxns(i);
    mytbl(i,2) = {flux(strcmp(model.rxns,myrxns{i}))};
    mytbl(i,3) = printRxnFormula(model, myrxns{i});
end
mytbl(:,4) = mat2cell(myfluxP,ones(length(myfluxP),1),1);
mytbl = sortrows(mytbl,4);
end