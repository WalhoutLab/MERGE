addpath /Users/xuhangli/Desktop/Walhout_Lab/conjoined_model_project/new_matlab_functions
%initCobraToolbox;
model = iCEL2model_xl('iCEL1314.xlsx');
%%
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
% add the dynamic constraints
model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('EXC9999',model.rxns)) = 1; 
model.S(end, strcmp('EXC0050',model.rxns)) = 0.01; 
model.csense(end+1) = 'L';
model.b(end+1) = 0;
model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('EXC9998',model.rxns)) = 1; 
model.S(end, strcmp('EXC0050',model.rxns)) = 0.01; 
model.csense(end+1) = 'L';
model.b(end+1) = 0;

model = changeObjective(model,'BIO0010');
optimizeCbModel(model,'max')
model = changeObjective(model,'RCC0005');
optimizeCbModel(model,'max')
%%
save('iCEL1314.mat','model');