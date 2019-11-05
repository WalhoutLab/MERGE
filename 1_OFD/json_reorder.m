%% to process json file from python environment
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

%%
myflux = OFD_all.Neurons;
%convert to cell
mynames = fieldnames(myflux);
for j = 1:length(mynames)
    flux_tmp(j) = myflux.(mynames{j});
end
[A B] = ismember(model.rxns,mynames);
myflux_reorder = flux_tmp(B)';
