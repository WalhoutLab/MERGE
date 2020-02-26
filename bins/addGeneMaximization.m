function [MILPproblem] = addGeneMaximization(MILPproblem, HGeneList, HRxnList,LRxnList, model)
%imply Vi_gene <= Vi_rxn1 + Vi_rxn2 + ... + Vi_rxnn, while Vi_gene is a
%binary variable
%the S matrix should be at left top corner
%HRxnList ... the high expressed reaction used in ordinary iMAT
%HGeneList ... the high expressed gene list for maximization

%% find the HRxn orders in input matrix 
HRindex0 = find(ismember(model.rxns, HRxnList));
%create A matrix
A = sparse(size(MILPproblem.A,1)+length(HGeneList),size(MILPproblem.A,2)+length(HGeneList));
[m,n,s] = find(MILPproblem.A);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end

for i = 1:length(HGeneList)
    myRxns = model.rxns(any(model.rxnGeneMat(:,strcmp(HGeneList(i),model.genes)),2),1);
    HRNames = HRxnList(ismember(HRxnList,myRxns));
    HRindex1 = find(ismember(model.rxns, HRNames));
    HRindex = find(ismember(HRindex0,HRindex1));
    A(size(MILPproblem.A,1)+i,length(model.rxns)+HRindex) = 1;
    A(size(MILPproblem.A,1)+i,length(model.rxns)+length(HRxnList)+length(LRxnList)+HRindex) = 1;
    A(size(MILPproblem.A,1)+i,size(MILPproblem.A,2)+i) = -1;
end

%create other inputs
lb = [MILPproblem.lb;zeros(length(HGeneList),1)];
ub = [MILPproblem.ub;ones(length(HGeneList),1)];
b = [MILPproblem.b;zeros(length(HGeneList),1)];
csense1(1:length(HGeneList)) = 'G';
csense = [MILPproblem.csense,csense1];
vartype1(1:length(HGeneList),1) = 'B';
vartype = [MILPproblem.vartype;vartype1];

%generate function output
MILPproblem.A = A;
MILPproblem.b = b;
MILPproblem.lb = lb;
MILPproblem.ub = ub;
MILPproblem.csense = csense;
MILPproblem.vartype = vartype;
end