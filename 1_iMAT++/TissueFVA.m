%% This is the main program to generate tissue FVA in TableS5
%% Load model
addpath ~/cobratoolbox/
addpath /share/pkg/gurobi/810/linux64/matlab/
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
fprintf('C. elegans Tissue Model FVA Started... \n');
fprintf('loading model and input parameters... \n');
totalTime = tic();
load('Tissue.mat');
% The loaded model is already constrained with default constraints and
% side/storage parameters. To modify a reaction constraint (RCC0005 as an example reaction):
% worm = changeRxnBounds(worm,'RCC0005',10,'l');
% To change the side nutrient proportion:
% worm.S(end-1, strcmp('EXC0050',worm.rxns)) = 0.01*9.6641833;%storage
% worm.S(end, strcmp('EXC0050',worm.rxns)) = 0.01*9.6641833; %side

% Reset some constraints to make the model ready for integration 
% Release the input constraints for integration 
model = changeRxnBounds(model,'EXC0050_L',-1000,'l');
model = changeRxnBounds(model,'EX00001_E',-1000,'l');
model = changeRxnBounds(model,'EX00007_E',-1000,'l');
model = changeRxnBounds(model,'EX00009_E',-1000,'l');
model = changeRxnBounds(model,'EX00080_E',-1000,'l');

model = changeRxnBounds(model,'RM00112_I',-1000,'l');
model = changeRxnBounds(model,'RM00112_X',-1000,'l');
model = changeRxnBounds(model,'RM00112_I',1000,'u');
model = changeRxnBounds(model,'RM00112_X',1000,'u');

model = changeRxnBounds(model,'DMN0033_I',0,'l');
model = changeRxnBounds(model,'DMN0033_X',0,'l');
model = changeRxnBounds(model,'DMN0033_I',1,'u');
model = changeRxnBounds(model,'DMN0033_X',1,'u');

model = changeRxnBounds(model,'RMC0005_I',0,'l');
model = changeRxnBounds(model,'RMC0005_X',0,'l');
model = changeRxnBounds(model,'RMC0005_I',0,'u');
model = changeRxnBounds(model,'RMC0005_X',0,'u');
%
% parseGPR takes hugh amount of time, so preparse and save result in model
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

%% Load flux thresholds (epsilons)
% The epsilon values are also supplied in a jason file. We convert it to a vector
% variable that follows the same order as model.rxns.
fname = 'epsilon.json';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
str = regexprep(str,'u"','"');
epsilon = jsondecode(str);
fields = fieldnames(epsilon);
epsilon_f = 0.01*ones(length(model.rxns),1); % the default epsilon is assumed as 0.01 in this study.
epsilon_r = 0.01*ones(length(model.rxns),1);
for i = 1: length(fields)
    if epsilon.(fields{i})(2) ~= 0
        epsilon_f(strcmp(model.rxns,fields{i})) = epsilon.(fields{i})(2);
    end
    if epsilon.(fields{i})(1) ~= 0
        epsilon_r(strcmp(model.rxns,fields{i})) = epsilon.(fields{i})(1);
    end
end

%% Define the par cluster
myCluster = parcluster('local');
myCluster.NumWorkers = 128;
saveProfile(myCluster);
parpool(20,'SpmdEnabled',false); % we tested this program in a 20-core server

%% Do FVA for the six non-intestinal tissues.
names = {'Body_wall_muscle','Glia','Gonad','Hypodermis','Neurons','Pharynx'};
Xrxns = model.rxns; % we calculate the FVA of both X tissue and I tissue
for i = 1:length(names)
    fprintf('now starting to calculate for %s... \n',names{i});
    fitTime = tic();
    targetRxns = Xrxns;
    parforFlag = 1;
    myFVA = struct(); %my context specific model
    load(['output/',names{i},'_speedmode_1.mat']);
    [myFVA.lb, myFVA.ub] = FVA_MILP(myCSM.MILP, model, targetRxns,parforFlag);
    save(['output/FVA/',names{i},'_speedmode_1.mat'],'myFVA');
    toc(fitTime);
end

%% Now do the intestine integration
% first get the average flux distribution of all six non-intestinal
% tissues.
for i = 1:length(names)
    load(['output/',names{i},'.mat']);
    eval([names{i},' = myCSM;']);
end
weights = [0.16481752424244808;0.19668335945561571;0.029394841229723516;0.055442983787911716;0.35973476877630634;0.19392652250799472];%needs to order by tissuesname
fluxSum_OFD = zeros(length(model.rxns),1);
for i = 1:length(names)
    eval(['fluxSum_OFD = ',names{i},'.OFD .* weights(i) + fluxSum_OFD;']);
end
% Collapse the X compartment in the dual model to save computational power
IntestineModel_OFD = collapseX(model,'X',fluxSum_OFD);
% Update parsed GPR
parsedGPR = GPRparser_xl(IntestineModel_OFD);% Extracting GPR data from model
IntestineModel_OFD.parsedGPR = parsedGPR;
% Update the epsilon list
epsilon_new_f = 0.01 * ones(length(IntestineModel_OFD.rxns),1);
epsilon_new_r = epsilon_new_f;
for i = 1: length(IntestineModel_OFD.rxns)
    if any(ismember(IntestineModel_OFD.rxns(i), model.rxns))
        epsilon_new_f(i) = epsilon_f(ismember(model.rxns,IntestineModel_OFD.rxns(i)));
        epsilon_new_r(i) = epsilon_r(ismember(model.rxns,IntestineModel_OFD.rxns(i)));
    end
end

load('output/Intestine_speedmode_1.mat');

% Set free all X exchange
myCSM.MILP.ub(cellfun(@(x) ~isempty(regexp(x, '_X$','once')) ,IntestineModel_OFD.rxns)) = 1000;
myCSM.MILP.lb(cellfun(@(x) ~isempty(regexp(x, '_X$','once')) ,IntestineModel_OFD.rxns)) = 0;

% Integration 
fprintf('now starting to calculate for %s... \n','Intestine');
fitTime = tic();
parforFlag=1;
% we only calculate the FVA of the I tissue
targetRxns = IntestineModel_OFD.rxns(cellfun(@(x) isempty(regexp(x, '_X$','once')) ,IntestineModel_OFD.rxns));
myFVA = struct(); %my context specific model
% OFD fitting
[myFVA.lb, myFVA.ub] = FVA_MILP(myCSM.MILP, IntestineModel_OFD, targetRxns,parforFlag);
save('output/FVA/Intestine_speedmode_1.mat','myFVA');
toc(fitTime);
fprintf('C. elegans Tissue Model FVA Completed! \n');
toc(totalTime);

%% For generating the list of reactions to block, see the step #3 of
% walkthrough_large_scale_FVA.m
