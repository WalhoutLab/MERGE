function []=DGD(i) %i gene in E.coli to all the worm gene
i = str2num(i);
addpath ~/cobratoolbox/
addpath /share/pkg/gurobi/810/linux64/matlab/
addpath ~/new_matlab_functions/
initCobraToolbox;

%% load model
load worm_new.mat
model = worm_new;
%% synthetic lethal
wormGenes = regexprep(model.genes,'_deleted','');
myWormGrowth_BIO0100 = sparse(1,length(wormGenes));
myWormGrowth_BIO0103 = sparse(1,length(wormGenes));
qryGene = wormGenes(i);
for j = 1:length(wormGenes)
    [model_del, ~, ~, ~] = deleteModelGenes(model, [qryGene;wormGenes(j)]);
    model_del = changeObjective(model_del,'BIO0100');
    opti = optimizeCbModel(model_del);
    if opti.stat == 1 %feasible
        myWormGrowth_BIO0100(j) = opti.obj;
    else
        myWormGrowth_BIO0100(j) = nan;
    end
    model_del = changeObjective(model_del,'BIO0103');
    opti = optimizeCbModel(model_del);
    if opti.stat == 1 %feasible
        myWormGrowth_BIO0103(j) = opti.obj;
    else
        myWormGrowth_BIO0103(j) = nan;
    end
j
end
save(['BIO0100Worm_',num2str(i),'.mat'],'myWormGrowth_BIO0100');
save(['BIO0103Worm_',num2str(i),'.mat'],'myWormGrowth_BIO0103');
end
