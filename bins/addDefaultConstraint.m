function [model] = addDefaultConstraint(model,type)
%%
if strcmp(type,'default')
    model = changeRxnBounds(model,'DMN0033',0,'l');
    model = changeRxnBounds(model,'DMN0033',1,'u');
    model = changeRxnBounds(model,'EX00001',-1000,'l');
    model = changeRxnBounds(model,'EX00007',-1000,'l');
    model = changeRxnBounds(model,'EX00009',-1000,'l');
    model = changeRxnBounds(model,'EX00080',-1000,'l');
    model = changeRxnBounds(model,'RCC0005',1,'l');
    model = changeRxnBounds(model,'RCC0005',1000,'u');
    model = changeRxnBounds(model,'RM00112',0,'l');
    model = changeRxnBounds(model,'RM00112',0,'u');
    model = changeRxnBounds(model,'RM04432',-1000,'l');
    model = changeRxnBounds(model,'RM04432',1000,'u');
    model = changeRxnBounds(model,'RMC0005',0,'u');
    model = changeRxnBounds(model,'RMC0005',0,'l');
    model = changeRxnBounds(model,'EXC0050',-1,'l');
    model = changeRxnBounds(model,'EXC0050',0,'u');
    %model = changeRxnBounds(model,'EX00187',-1000,'l');
    model = changeObjective(model,'BIO0010');
elseif strcmp(type,'nutritionalFree@1')
    model.lb(cellfun(@(x) strcmp(x(1:2),'EX'),model.rxns)) = -1;
    model = changeObjective(model,'BIO0010');
    model = changeRxnBounds(model,'DMN0033',0,'l');
    model = changeRxnBounds(model,'DMN0033',0,'u');
    model = changeRxnBounds(model,'EX00001',-1000,'l');
    model = changeRxnBounds(model,'EX00007',-1000,'l');
    model = changeRxnBounds(model,'EX00009',-1000,'l');
    model = changeRxnBounds(model,'EX00080',-1000,'l');
    model = changeRxnBounds(model,'RCC0005',1,'l');
    model = changeRxnBounds(model,'RCC0005',1000,'u');
    model = changeRxnBounds(model,'RM00112',0,'l');
    model = changeRxnBounds(model,'RM00112',0,'u');
    model = changeRxnBounds(model,'RM04432',-1000,'l');
    model = changeRxnBounds(model,'RM04432',1000,'u');
    model = changeRxnBounds(model,'RMC0005',0,'u');
    model = changeRxnBounds(model,'RMC0005',0,'l');
    model = changeRxnBounds(model,'EXC0050',-1,'l');
    model = changeRxnBounds(model,'EXC0050',0,'u');
    model = changeRxnBounds(model,'EX00187',-1000,'l');
    model = changeObjective(model,'BIO0010');
elseif strcmp(type,'nutritionalFree@1000')
    model.lb(cellfun(@(x) strcmp(x(1:2),'EX'),model.rxns)) = -1000;
    model = changeObjective(model,'BIO0010');
    model = changeRxnBounds(model,'DMN0033',0,'l');
    model = changeRxnBounds(model,'DMN0033',1000,'u'); % be careful about this. could be 1
    model = changeRxnBounds(model,'EX00001',-1000,'l');
    model = changeRxnBounds(model,'EX00007',-1000,'l');
    model = changeRxnBounds(model,'EX00009',-1000,'l');
    model = changeRxnBounds(model,'EX00080',-1000,'l');
    model = changeRxnBounds(model,'RCC0005',1,'l');
    model = changeRxnBounds(model,'RCC0005',1000,'u');
    model = changeRxnBounds(model,'RM00112',-1000,'l');
    model = changeRxnBounds(model,'RM00112',1000,'u');
    model = changeRxnBounds(model,'RM04432',-1000,'l');
    model = changeRxnBounds(model,'RM04432',1000,'u');
    model = changeRxnBounds(model,'RMC0005',0,'u');
    model = changeRxnBounds(model,'RMC0005',0,'l');
    model = changeRxnBounds(model,'EXC0050',-1000,'l');
    model = changeRxnBounds(model,'EXC0050',0,'u');
    model = changeRxnBounds(model,'EX00187',-1000,'l');
    model = changeRxnBounds(model,'EXC0216',0,'l'); % special molecular that can only be produced
    model = changeRxnBounds(model,'EX01060',0,'l'); % special molecular that can only be produced
    model = changeRxnBounds(model,'EX18125',0,'l'); % special molecular that can only be produced
    model = changeRxnBounds(model,'EX18126',0,'l'); % special molecular that can only be produced

elseif strcmp(type,'default-Conjioned') 
ddd    model = changeRxnBounds(model,'DMN0033_W',0,'l');
    model = changeRxnBounds(model,'DMN0033_W',0,'u');
    model = changeRxnBounds(model,'EX00001_W',-1000,'l');
    model = changeRxnBounds(model,'EX00007_W',-1000,'l');
    model = changeRxnBounds(model,'EX00009_W',-1000,'l');
    model = changeRxnBounds(model,'EX00080_W',-1000,'l');
    model = changeRxnBounds(model,'RCC0005_W',1,'l');
    model = changeRxnBounds(model,'RCC0005_W',1000,'u');
    model = changeRxnBounds(model,'RM00112_W',0,'l');
    model = changeRxnBounds(model,'RM00112_W',0,'u');
    model = changeRxnBounds(model,'RM04432_W',-1000,'l');
    model = changeRxnBounds(model,'RM04432_W',1000,'u');
    model = changeRxnBounds(model,'RMC0005_W',0,'u');
    model = changeRxnBounds(model,'RMC0005_W',0,'l');
    model = changeRxnBounds(model,'EX00187_W',-1000,'l');
    % constrians for e.coli is not added!
    model = changeObjective(model,'BIO0010_W');
elseif strcmp(type,'yeast-SL') 
    model = deleteModelGenes(model,'YEL063C');
    model = deleteModelGenes(model,'YNL268W');
    model = deleteModelGenes(model,'YEL021W');
    model = deleteModelGenes(model,'YCL018W');
    model = deleteModelGenes(model,'YLR303W');
% five genotype for the yeast in brendas data
% close all exchange and only set the defined ones
    model = changeRxnBounds(model,model.rxns(findExcRxns(model)),0,'l');
    %now define the media
    model = changeRxnBounds(model,'r_1604',-0.000002,'l');
    model = changeRxnBounds(model,'r_1639',-3.01,'l');
    model = changeRxnBounds(model,'r_1873',-0.36,'l');
    model = changeRxnBounds(model,'r_1880',-0.36,'l');
    model = changeRxnBounds(model,'r_1881',-0.36,'l');
    model = changeRxnBounds(model,'r_1671',-0.00000142,'l');
    model = changeRxnBounds(model,'r_1883',-0.36,'l');
    model = changeRxnBounds(model,'r_1861',-1000,'l');
    model = changeRxnBounds(model,'r_1714',-22.6,'l');
    model = changeRxnBounds(model,'r_1891',-0.36,'l');
    model = changeRxnBounds(model,'r_1889',-3.6,'l');
    model = changeRxnBounds(model,'r_1810',-0.36,'l');
    model = changeRxnBounds(model,'r_1897',-0.36,'l');
    model = changeRxnBounds(model,'r_1947',-0.11,'l');
    model = changeRxnBounds(model,'r_2020',-4.44,'l');
    model = changeRxnBounds(model,'r_1899',-1.8,'l');
    model = changeRxnBounds(model,'r_1902',-0.36,'l');
    model = changeRxnBounds(model,'r_2049',-0.75,'l');
    model = changeRxnBounds(model,'r_1967',-0.000002,'l');
    model = changeRxnBounds(model,'r_1992',-6.3,'l');
    model = changeRxnBounds(model,'r_1903',-0.36,'l');
    model = changeRxnBounds(model,'r_2005',-0.89,'l');
    model = changeRxnBounds(model,'r_1548',-0.0002,'l');
    model = changeRxnBounds(model,'r_1904',-0.36,'l');
    model = changeRxnBounds(model,'r_2038',-0.00092,'l');
    model = changeRxnBounds(model,'r_1906',-0.36,'l');
    model = changeRxnBounds(model,'r_2060',-100,'l');
    model = changeRxnBounds(model,'r_2067',-0.0032,'l');
    model = changeRxnBounds(model,'r_1911',-0.36,'l');
    model = changeRxnBounds(model,'r_1912',-0.36,'l');
    model = changeRxnBounds(model,'r_1913',-0.36,'l');
    model = changeRxnBounds(model,'r_2090',-3.63,'l');
    model = changeRxnBounds(model,'r_1914',-0.36,'l');
    model = changeRxnBounds(model,'r_2100',-1000,'l');%h2o
    model = changeRxnBounds(model,'r_1832',-1000,'l');%h
    model = changeRxnBounds(model,'r_4593',-1000,'l');%cl
    model = changeRxnBounds(model,'r_4594',-1,'l');%cu
    model = changeRxnBounds(model,'r_4595',-1,'l');%mn
    model = changeRxnBounds(model,'r_4596',-1,'l');%zn
    model = changeRxnBounds(model,'r_4597',-1,'l');%mg
    model = changeRxnBounds(model,'r_4600',-1,'l');%ca
    %user defined according to BD difco
    model = changeRxnBounds(model,'r_1792',-0.000002,'l');%folic acid
    model = changeRxnBounds(model,'r_2028',-0.000002,'l');%folic acid

    % constrians for e.coli is not added!
    model = changeObjective(model,'r_2111');
end
end
