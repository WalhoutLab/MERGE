addpath /Users/xuhangli/Desktop/Walhout_Lab/conjoined_model_project/new_matlab_functions
%initCobraToolbox;
% the <- in rxns table should be change to <==> manually
% ub/lb header to "Lower bound	Upper bound"
% change met header
% add two psudo met to the met table
% the rxns table and met table should be merged into one excel table
model = iCEL2model_xl('./input/model.xlsx');
%% add the dynamic constraints
sideR=0.02;
sideR_ind=0.005;
storageR=0.01;
bacMW=966.28583751;
% dividedBy=100; this is in the stoicheomitry matrix
% add individual met control and scale by the factor
sideRxns =  model.rxns(model.S(strcmp(model.mets,'sideMet[e]'),:)~=0);
sideRxns(strcmp(sideRxns,'EXC9998_E')) = [];
for i = 1:length(sideRxns)
    model.S(end+1,:) = zeros(1,length(model.rxns));
    model.S(end, strcmp(sideRxns{i},model.rxns)) = model.S(strcmp(model.mets,'sideMet[e]'),strcmp(sideRxns{i},model.rxns)); 
    model.S(end, strcmp('EXC0050_L',model.rxns)) = sideR_ind*bacMW*0.01;%the 0.01 here is the 1/devidedBy 
    model.csense(end+1) = 'L';
    model.b(end+1) = 0;
    model.BiGG(end+1) = {'NA'};
    model.metCharges(end+1) = 0;
    model.metFormulas(end+1) = {'NA'};
    model.metNames(end+1) = {'Nonsense Metabolites which is a constraint'};
    model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};
end

% add overall side/storage control
model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('EXC9999_E',model.rxns)) = 1; 
model.S(end, strcmp('EXC0050_L',model.rxns)) = storageR*bacMW*0.01; 
model.csense(end+1) = 'L';
model.b(end+1) = 0;
model.BiGG(end+1) = {'NA'};
model.metCharges(end+1) = 0;
model.metFormulas(end+1) = {'NA'};
model.metNames(end+1) = {'Nonsense Metabolites which is a constraint'};
model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};

model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('EXC9998_E',model.rxns)) = 1; 
model.S(end, strcmp('EXC0050_L',model.rxns)) = sideR*bacMW*0.01; 
model.csense(end+1) = 'L';
model.b(end+1) = 0;
model.BiGG(end+1) = {'NA'};
model.metCharges(end+1) = 0;
model.metFormulas(end+1) = {'NA'};
model.metNames(end+1) = {'Nonsense Metabolites which is a constraint'};
model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};

model = changeObjective(model,'BIO0010_I');
opt = optimizeCbModel(model,'max')
model = changeObjective(model,'RCC0005_X');
optimizeCbModel(model,'max')
%% 
save('Tissue.mat','model');
%% flux minimization
test = changeRxnBounds(model,'BIO0010_I',opt.obj-1e-5,'l');
minflux = minimizeModelFlux_XL(test)

%% the following makes the intestine model
%initCobraToolbox;
% the <- in rxns table should be change to <==> manually
% ub/lb header to "Lower bound	Upper bound"
% change met header
% add two psudo met to the met table
% the rxns table and met table should be merged into one excel table
model = iCEL2model_xl('./input/intestine.xlsx');
%% add the dynamic constraints
sideR=0.02;
sideR_ind=0.005;
storageR=0.01;
bacMW=966.28583751;
% dividedBy=100; this is in the stoicheomitry matrix
% add individual met control and scale by the factor
sideRxns =  model.rxns(model.S(strcmp(model.mets,'sideMet[e]'),:)~=0);
sideRxns(strcmp(sideRxns,'EXC9998_E')) = [];
for i = 1:length(sideRxns)
    model.S(end+1,:) = zeros(1,length(model.rxns));
    model.S(end, strcmp(sideRxns{i},model.rxns)) = model.S(strcmp(model.mets,'sideMet[e]'),strcmp(sideRxns{i},model.rxns)); 
    model.S(end, strcmp('EXC0050_L',model.rxns)) = sideR_ind*bacMW*0.01;%the 0.01 here is the 1/devidedBy 
    model.csense(end+1) = 'L';
    model.b(end+1) = 0;
    model.BiGG(end+1) = {'NA'};
    model.metCharges(end+1) = 0;
    model.metFormulas(end+1) = {'NA'};
    model.metNames(end+1) = {'Nonsense Metabolites which is a constraint'};
    model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};
end

% add overall side/storage control
model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('EXC9999_E',model.rxns)) = 1; 
model.S(end, strcmp('EXC0050_L',model.rxns)) = storageR*bacMW*0.01; 
model.csense(end+1) = 'L';
model.b(end+1) = 0;
model.BiGG(end+1) = {'NA'};
model.metCharges(end+1) = 0;
model.metFormulas(end+1) = {'NA'};
model.metNames(end+1) = {'Nonsense Metabolites which is a constraint'};
model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};

model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('EXC9998_E',model.rxns)) = 1; 
model.S(end, strcmp('EXC0050_L',model.rxns)) = sideR*bacMW*0.01; 
model.csense(end+1) = 'L';
model.b(end+1) = 0;
model.BiGG(end+1) = {'NA'};
model.metCharges(end+1) = 0;
model.metFormulas(end+1) = {'NA'};
model.metNames(end+1) = {'Nonsense Metabolites which is a constraint'};
model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};

model = changeObjective(model,'BIO0010_I');
opt = optimizeCbModel(model,'max')
model = changeObjective(model,'RCC0005_X');
optimizeCbModel(model,'max')
%% 
save('Intestine.mat','model');
