%% to compare safak's flux distribution with matlab result
safak = readtable('input/pdf_XuhangsParameters.txt');
%%
names = fieldnames(fittingQuality.D_Rlow)
mysum = 0
for i = 1:length(names)
    mysum = mysum + abs(fittingQuality.D_Rlow.(names{i}));
end
mysum
%%
sum(abs(safak.FLUX))
sum(abs(PFD))
sum(abs(safak.FLUX(ismember(model.rxns,names))))
sum(abs(PFD(ismember(model.rxns,names))))