function [rxns_annotated, subsystemInfo] = annotateRxnSet(rxns,model)
% annotate the subsystems and calculate the p-value of enrichment 
%rxn        formula     subsystem
%subsystem  p-value
rxns_annotated(:,1) = rxns;
for i = 1:length(rxns)
    rxns_annotated(i,2) = printRxnFormula(model,'rxnAbbrList',rxns(i),'printFlag',0);
    rxns_annotated(i,3) = model.subSystems(strcmp(model.rxns,rxns{i}));
end
% calculate p-value
subsystemInfo(:,1) = unique(model.subSystems(ismember(model.rxns,rxns)));
for i = 1:length(subsystemInfo(:,1))
    mysys = subsystemInfo(i,1);
    n_obs = sum(ismember(model.subSystems(ismember(model.rxns,rxns)),mysys));
    n_sample = length(rxns);
    n_pop = length(model.rxns);
    n_trueInPop = sum(ismember(model.subSystems,mysys));
    subsystemInfo{i,2} = 1-hygecdf(n_obs-1,n_pop,n_trueInPop,n_sample);
end
subsystemInfo = sortrows(subsystemInfo,2);
end