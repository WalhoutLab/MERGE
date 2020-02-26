fname = './input/OFD_all.json';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
str = regexprep(str,'u"','"');
str = regexprep(str,'set(','');
str = regexprep(str,')','');
OFD_all = jsondecode(str);

names = fieldnames(tissues);
names = names(1:end-1);
weights = [0.179980982;0.086522532;0.036874975;0.138793829;0.372097477000000;0.185730205];%needs to order by tissuesname
%weights = [0.179980982;0.086522532;0.036874975;0.185730205;0.372097477000000;0.138793829];
fluxSum = zeros(length(model.rxns),1);
flux_tmp = [];
for i = 1:length(names)
    myflux = OFD_all.(names{i});
    %convert to cell
    mynames = fieldnames(myflux);
    for j = 1:length(mynames)
        flux_tmp(j) = myflux.(mynames{j});
    end
    [A B] = ismember(model.rxns,mynames);
    myflux_reorder = flux_tmp(B)';
    fluxSum = myflux_reorder .* weights(i) + fluxSum;
end
%%
IntestineModel = collapseX(model,'X',fluxSum);