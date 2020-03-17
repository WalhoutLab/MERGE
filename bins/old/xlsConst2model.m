function model_constrained  = xlsConst2model(fileName, model)
%the input xls should have 3 columns as name, lower and upper bound 
% initialization with default values
[num,txt,all] = xlsread(fileName, 'Constrain');
upper = num(:,2);
lower = num(:,1);
names1 = all(2:end,1);
model_constrained = changeRxnBounds(model, names1, upper, 'u');
model_constrained = changeRxnBounds(model_constrained, names1, lower, 'l');
end