names = {'Hypodermis','Gonad','Glia','Pharynx','Body_wall_muscle','Neurons','Intestine'};
ref = readtable('./../../../ofd_table.tsv','FileType','text');
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
%% FVA
load('output/FVA/Hypodermis.mat');
load('output/Hypodermis.mat');
rxnID =  model.rxns(cellfun(@(x) ~isempty(regexp(x,'_X$','once')),model.rxns));
levels_f = zeros(length(rxnID),1);
levels_r = zeros(length(rxnID),1);
for i = 1:length(rxnID)
    myrxn = rxnID{i};
    if myCSM.OFD(strcmp(model.rxns,myrxn)) > 1e-5
        levels_f(i) = 1;
    elseif myCSM.OFD(strcmp(model.rxns,myrxn)) > 1e-7 && myFVA.lb(i) > 1e-7
        levels_f(i) = 1;
    elseif myFVA.ub(i) > max(epsilon_f(strcmp(model.rxns,myrxn))-1e-5,1e-5)
        levels_f(i) = -1;
    else
        if myFVA.ub(i) == 0
            levels_f(i) = -3;
        else
            levels_f(i) = -2;
        end
    end
    
    if -myCSM.OFD(strcmp(model.rxns,myrxn)) > 1e-5
        levels_r(i) = 1;
    elseif -myCSM.OFD(strcmp(model.rxns,myrxn)) > 1e-7 && -myFVA.ub(i) > 1e-7
        levels_r(i) = 1;
    elseif -myFVA.lb(i) > max(epsilon_r(strcmp(model.rxns,myrxn))-1e-5,1e-5)
        levels_r(i) = -1;
    else
        if myFVA.lb(i) == 0
            levels_r(i) = -3;
        else
            levels_r(i) = -2;
        end
    end
end
levels_f(levels_f==-2)=-3;
levels_r(levels_r==-2)=-3;
%%      
confidenceTbl = readtable('confidenceTable.tsv','FileType','text','ReadRowNames',true);
confidenceTbl.Properties.RowNames = cellfun(@(x) [x(1:end-1),'_',x(end)],confidenceTbl.Properties.RowNames,'UniformOutput',false);
totalN = 0;
matchN = 0;
confidenceTbl.Hypodermis(confidenceTbl.Hypodermis==-2) = -3;
for i = 1:length(rxnID)
    if any(strcmp(confidenceTbl.Properties.RowNames,[rxnID{i},'_f']))
        totalN = totalN+1;
        if confidenceTbl.Hypodermis(strcmp(confidenceTbl.Properties.RowNames,[rxnID{i},'_f'])) == levels_f(i)
            matchN= matchN+1;
        else
            break
        end
    end
    if any(strcmp(confidenceTbl.Properties.RowNames,[rxnID{i},'_r']))
        totalN = totalN+1;
        if confidenceTbl.Hypodermis(strcmp(confidenceTbl.Properties.RowNames,[rxnID{i},'_r'])) == levels_r(i)
            matchN= matchN+1;
        end
    end
end
    
    