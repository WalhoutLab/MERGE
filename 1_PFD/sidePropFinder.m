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
% note side rxn like UP18125 ==> may represent a bulk content!
wasteDW = [];
Nfit = [];
bacUp = [];
for sideP = 0.04:0.01:0.25
    [FluxDistribution_daf2,N_highFit_daf2,N_zeroFit_daf2,minLow_daf2,minTotal_daf2,~,wasteDW(end+1),Hgenes,Lrxns] = autoIntegration_latent(worm,true,0.01,sideP,epsilon_f,epsilon_r,capacity_f,capacity_r, 10, ExpCatag);
    Nfit(end+1) = N_highFit_daf2+N_zeroFit_daf2;
    bacUp(end+1) = FluxDistribution_daf2(strcmp(worm.rxns,'EXC0050'));
    sideP
    out = fopen('./output.txt','a+');%write the infeasible cases to a log file
    fprintf(out,'%f\t%.2f\t%d\t%.2f\n',sideP,wasteDW(end),Nfit(end),bacUp(end));
    fclose(out);
end
%FluxDistribution(strcmp(worm.rxns,'EXC0051'))
