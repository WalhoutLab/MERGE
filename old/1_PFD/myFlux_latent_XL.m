%% load model
addpath ~/cobratoolbox/
addpath /share/pkg/gurobi/810/linux64/matlab/
addpath ~/new_matlab_functions/
initCobraToolbox;
load('iCEL1314.mat');
worm = addDefaultConstraint(model,'default');
worm = changeRxnBounds(worm,'RCC0005',10,'l');
% change the side proportion
worm.S(end-1, strcmp('EXC0050',worm.rxns)) = 0.01*9.6641833;
worm.S(end, strcmp('EXC0050',worm.rxns)) = 0.01*9.6641833;
% make the epsilon vector
[epsilon_f,epsilon_r,capacity_f,capacity_r] = makeEpsilonSeq(worm,worm.rxns,0.01,0.5);

%% daf2
load refined_daf2.mat
% map to iCEL names
fname = 'geneNameTable.txt';
str = fileread(fname);
for i = 1:length(str)
    if  str(i) == "'"
        str(i) = '"';
    end
end
GeneNameTable = jsondecode(str);
for i = 1:length(ExpCatag.dynamic)
    ExpCatag.dynamic{i} = GeneNameTable.(ExpCatag.dynamic{i});
end
for i = 1:length(ExpCatag.high)
    ExpCatag.high{i} = GeneNameTable.(ExpCatag.high{i});
end
for i = 1:length(ExpCatag.low)
    ExpCatag.low{i} = GeneNameTable.(ExpCatag.low{i});
end
for i = 1:length(ExpCatag.zero)
    ExpCatag.zero{i} = GeneNameTable.(ExpCatag.zero{i});
end

%% filter and write flux
myFlux = struct();
[myFlux.OFD,myFlux.N_highFit,myFlux.N_zeroFit,myFlux.minLow,myFlux.minTotal,myFlux.OpenGene,myFlux.wasteDW,myFlux.HGenes,myFlux.RLNames,myFlux.latentRxns,myFlux.sensitivity,myFlux.PFD] = autoIntegration_latent(worm,true,0.01,0.05,epsilon_f,epsilon_r,capacity_f,capacity_r, 24, ExpCatag);
myFlux.OFD = fix(myFlux.OFD .* 1e7) ./ 1e7;
myFlux.PFD = fix(myFlux.PFD .* 1e7) ./ 1e7;
save('daf2Flux.mat','myFlux');
%% wt
load refined_N2.mat
% map to iCEL names
fname = 'geneNameTable.txt';
str = fileread(fname);
for i = 1:length(str)
    if  str(i) == "'"
        str(i) = '"';
    end
end
GeneNameTable = jsondecode(str);
for i = 1:length(ExpCatag.dynamic)
    ExpCatag.dynamic{i} = GeneNameTable.(ExpCatag.dynamic{i});
end
for i = 1:length(ExpCatag.high)
    ExpCatag.high{i} = GeneNameTable.(ExpCatag.high{i});
end
for i = 1:length(ExpCatag.low)
    ExpCatag.low{i} = GeneNameTable.(ExpCatag.low{i});
end
for i = 1:length(ExpCatag.zero)
    ExpCatag.zero{i} = GeneNameTable.(ExpCatag.zero{i});
end
%%
myFlux = struct();
[myFlux.OFD,myFlux.N_highFit,myFlux.N_zeroFit,myFlux.minLow,myFlux.minTotal,myFlux.OpenGene,myFlux.wasteDW,myFlux.HGenes,myFlux.RLNames,myFlux.latentRxns,myFlux.sensitivity,myFlux.PFD] = autoIntegration_latent(worm,true,0.01,0.05,epsilon_f,epsilon_r,capacity_f,capacity_r, 8, ExpCatag);
myFlux.OFD = fix(myFlux.OFD .* 1e7) ./ 1e7;
myFlux.PFD = fix(myFlux.PFD .* 1e7) ./ 1e7;
save('N2Flux.mat','myFlux');
%% daf16
% load refined_daf16.mat
% % map to iCEL names
% fname = 'geneNameTable.txt';
% str = fileread(fname);
% for i = 1:length(str)
%     if  str(i) == "'"
%         str(i) = '"';
%     end
% end
% GeneNameTable = jsondecode(str);
% for i = 1:length(ExpCatag.dynamic)
%     ExpCatag.dynamic{i} = GeneNameTable.(ExpCatag.dynamic{i});
% end
% for i = 1:length(ExpCatag.high)
%     ExpCatag.high{i} = GeneNameTable.(ExpCatag.high{i});
% end
% for i = 1:length(ExpCatag.low)
%     ExpCatag.low{i} = GeneNameTable.(ExpCatag.low{i});
% end
% for i = 1:length(ExpCatag.zero)
%     ExpCatag.zero{i} = GeneNameTable.(ExpCatag.zero{i});
% end
% %% calculate the flux distribution
% myFlux = struct();
% [myFlux.OFD,myFlux.N_highFit,myFlux.N_zeroFit,myFlux.minLow,myFlux.minTotal,myFlux.OpenGene,myFlux.wasteDW,myFlux.HGenes,myFlux.RLNames,myFlux.latentRxns,myFlux.sensitivity,myFlux.PFD] = autoIntegration_latent(worm,true,0.01,sideP,epsilon_f,epsilon_r,capacity_f,capacity_r, 10, ExpCatag);
% myFlux.OFD = fix(myFlux.OFD .* 1e7) ./ 1e7;
% myFlux.PFD = fix(myFlux.PFD .* 1e7) ./ 1e7;
% save('daf16Flux.mat','daf16Flux');
%% double
load refined_double.mat
% map to iCEL names
fname = 'geneNameTable.txt';
str = fileread(fname);
for i = 1:length(str)
    if  str(i) == "'"
        str(i) = '"';
    end
end
GeneNameTable = jsondecode(str);
for i = 1:length(ExpCatag.dynamic)
    ExpCatag.dynamic{i} = GeneNameTable.(ExpCatag.dynamic{i});
end
for i = 1:length(ExpCatag.high)
    ExpCatag.high{i} = GeneNameTable.(ExpCatag.high{i});
end
for i = 1:length(ExpCatag.low)
    ExpCatag.low{i} = GeneNameTable.(ExpCatag.low{i});
end
for i = 1:length(ExpCatag.zero)
    ExpCatag.zero{i} = GeneNameTable.(ExpCatag.zero{i});
end

%% calculate the flux distribution
myFlux = struct();
[myFlux.OFD,myFlux.N_highFit,myFlux.N_zeroFit,myFlux.minLow,myFlux.minTotal,myFlux.OpenGene,myFlux.wasteDW,myFlux.HGenes,myFlux.RLNames,myFlux.latentRxns,myFlux.sensitivity,myFlux.PFD] = autoIntegration_latent(worm,true,0.01,0.05,epsilon_f,epsilon_r,capacity_f,capacity_r, 8, ExpCatag);
myFlux.OFD = fix(myFlux.OFD .* 1e7) ./ 1e7;
myFlux.PFD = fix(myFlux.PFD .* 1e7) ./ 1e7;
save('doubleFlux.mat','myFlux');
% %% save the model parameter
% constraints = constraints2JSON(worm);
% fid=fopen('Constraints_XL.json','w');
% fprintf(fid, constraints);
% fclose(fid);
% 
% epsilons = epsilon2JSON(worm,epsilon_r,epsilon_f);
% fid=fopen('epsilons_XL.json','w');
% fprintf(fid, epsilons);
% fclose(fid);
% 
% genesets = geneSet2JSON(ExpCatag);
% fid=fopen('genesets_XL.json','w');
% fprintf(fid, genesets);
% fclose(fid);

