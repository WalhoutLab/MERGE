fname = './OFD.json';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
str = regexprep(str,'u"','"');
str = regexprep(str,'set(','');
str = regexprep(str,')','');
OFD_safak = jsondecode(str);
%%
myflux =OFD_safak.Body_wall_muscle;
%convert to cell
mynames = fieldnames(myflux);
flux_tmp = [];
for j = 1:length(mynames)
    flux_tmp(j) = myflux.(mynames{j});
end
mynames(strcmp(mynames,'EXEC9999_X')) =  {'EXC9999_X'};
mynames(strcmp(mynames,'EXEC9998_X')) = {'EXC9998_X'};
[A B] = ismember(model.rxns,mynames);
myflux_reorder = flux_tmp(B)';
%%
plot(myflux_reorder,myCSM.OFD,'.')
sum(abs(myflux_reorder))
sum(abs(myCSM.OFD))
%%
%OFD
names = fieldnames(tissues);
names = names(1:end-1);
tbl = table('Size',[length(model.rxns),6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',names','RowNames',model.rxns);
for i = 1:length(names)
    load(['./output/',names{i},'.mat']);
    tbl(:,i) = num2cell(myCSM.OFD);
end
writetable(tbl,'./output/allflux_OFD.csv','WriteRowNames', true);

%PFD
names = fieldnames(tissues);
names = names(1:end-1);
tbl = table('Size',[length(model.rxns),6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',names','RowNames',model.rxns);
for i = 1:length(names)
    load(['./output/',names{i},'.mat']);
    tbl(:,i) = num2cell(myCSM.PFD);
end
writetable(tbl,'./output/allflux_PFD.csv','WriteRowNames', true);

%intestine
names = {'OFD','PFD'};
tbl = table('Size',[length(IntestineModel.rxns),2],'VariableTypes',{'double','double'},'VariableNames',names,'RowNames',IntestineModel.rxns);
for i = 1:2
    load(['./output/Intestine.mat']);
    tbl(:,i) = num2cell(myCSM.(names{i}));
end
writetable(tbl,'./output/intestineFlux.csv','WriteRowNames', true);
%%
%%
%PFD
names = fieldnames(tissues);
names = names(1:end-1);
tbl = table('Size',[length(model.rxns),6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',names','RowNames',model.rxns);
for i = 1:length(names)
    load(['./output/',names{i},'.mat']);
    tbl(:,i) = num2cell(myCSM.PFD);
end
writetable(tbl,'./output/allflux_PFD.csv','WriteRowNames', true);

%intestine
names = {'PFD'};
tbl = table('Size',[length(IntestineModel.rxns),1],'VariableTypes',{'double'},'VariableNames',names,'RowNames',IntestineModel.rxns);
for i = 1
    load(['./output/Intestine.mat']);
    tbl(:,i) = num2cell(myCSM.(names{i}));
end
writetable(tbl,'./output/intestineFlux.csv','WriteRowNames', true);



