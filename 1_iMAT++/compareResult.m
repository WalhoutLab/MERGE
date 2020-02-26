names = {'Hypodermis','Gonad','Glia','Pharynx','Body_wall_muscle','Neurons','Intestine'};
ref = readtable('./../input/final/ofd_table.tsv','FileType','text');
for i = 1:length(names)-1
    load(['output/',names{i},'.mat']);
    myFlux = myCSM.OFD;
    refFlux = ref.(names{i});
    absDiff = sum(abs(myFlux - refFlux))
end
%% PFD
ref = readtable('./../input/pfd_table.tsv','FileType','text');
for i = 1:length(names)-1
    load(['output/',names{i},'.mat']);
    myFlux = myCSM.PFD;
    refFlux = ref.(names{i});
    absDiff = sum(abs(myFlux - refFlux));
    sum(abs(myFlux)) - sum(abs(refFlux))
end

%% show reaction 
myrxn = 'BIO0010';
ref = readtable('./../input/ofd_table.tsv','FileType','text');
for i = 1:length(names)-1
    load(['output/',names{i},'.mat']);
    myFlux = myCSM.OFD;
    refFlux = ref.(names{i});
    tbl(i,1) = myFlux(strcmp(model.rxns,[myrxn,'_X']));
    tbl(i,2) = refFlux(strcmp(model.rxns,[myrxn,'_X']));
end
load(['output/Intestine.mat']);
myFlux = myCSM.OFD;
refFlux = ref.Intestine;
tbl(i+1,1) = myFlux(strcmp(IntestineModel.rxns,[myrxn,'_I']));
tbl(i+1,2) = refFlux(strcmp(model.rxns,[myrxn,'_I']));
%% save all result to csv table
names = {'Hypodermis','Gonad','Glia','Pharynx','Body_wall_muscle','Neurons','Intestine'};
OFD = table();
PFD = table();
Intestine = table();
for i = 1:length(names)-1
    load(['output/',names{i},'.mat']);
    OFD.(names{i}) = myCSM.OFD;
    PFD.(names{i}) = myCSM.PFD;
end
OFD.Properties.RowNames = model.rxns;
PFD.Properties.RowNames = model.rxns;
writetable(OFD,'output/OFD.csv','Delimiter',',','WriteRowNames',true);
writetable(PFD,'output/PFD.csv','Delimiter',',','WriteRowNames',true);

load(['output/Intestine.mat']);
Intestine.PFD = myCSM.PFD;
Intestine.OFD = myCSM.OFD;
Intestine.Properties.RowNames = IntestineModel.rxns;
writetable(Intestine,'output/Intestine.csv','Delimiter',',','WriteRowNames',true);

