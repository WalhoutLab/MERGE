%% save all result to csv table
names = {'Hypodermis','Gonad','Glia','Pharynx','Body_wall_muscle','Neurons','Intestine'};
PFD = table();
Intestine = table();
for i = 1:length(names)-1
    load(['output/',names{i},'.mat']);
    PFD.(names{i}) = myCSM.PFD;
end
PFD.Properties.RowNames = model.rxns;
writetable(PFD,'output/PFD.csv','Delimiter',',','WriteRowNames',true);

load(['output/Intestine.mat']);
Intestine.PFD = myCSM.PFD;
Intestine.Properties.RowNames = IntestineModel.rxns;
writetable(Intestine,'output/Intestine.csv','Delimiter',',','WriteRowNames',true);

