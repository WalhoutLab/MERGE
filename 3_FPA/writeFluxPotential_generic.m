addpath ./model
addpath ./programs/2.refinement/
addpath ./programs/3.FluxFitting/
% prepare files
load('iCEL1315.mat');
model = addDefaultConstraint(model,'default');
% parseGPR takes hugh amount of time, so preparse and integrate with model
% here 
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

biomassList = readtable('biomassPrecursors.csv','ReadVariableNames',false,'Delimiter','\t');
biomassList = biomassList.Var1;

% expression files
fname = 'geneNameTable.txt';
str = fileread(fname);
for i = 1:length(str)
    if  str(i) == "'"
        str(i) = '"';
    end
end
GeneNameTable = jsondecode(str);
%distance = readtable('./safakFiles/iCEL1315_distance.txt');
load('./programs/1.initialCategory/TPM_FC.mat');
for i = 2:length(TPM_FC(:,1))
    TPM_FC{i,1} = GeneNameTable.(TPM_FC{i,1});
end
master_expression = {};
for i = 1:4
    expression = struct();
    expression.genes = TPM_FC(2:end,1);
    expression.value = cell2mat(TPM_FC(2:end,i+1));
    master_expression{i} = expression;
end
% load the distance matrix
load distanceMatrix.mat
% load the special penalties
manualPenalty = table2cell(readtable('manualPenalty.csv','ReadVariableNames',false,'Delimiter',','));
% load the special distance
manualDist = table2cell(readtable('manualDist.csv','ReadVariableNames',false,'Delimiter',','));

%% safak expression
sfTPM = table2cell(readtable('expression.tsv','FileType','text','Delimiter','\t'));
master_expression = {};
for i = 1:4
    expression = struct();
    expression.genes = sfTPM(:,1);
    expression.value = cell2mat(sfTPM(:,i+1));
    master_expression{i} = expression;
end

%% convert the file to a matrix
% distance = table2cell(distance);
% labels = unique(distance(:,1));
% distance_num = cell2mat(distance(:,3));
% distMat = nan(length(labels),length(labels));
% for i = 1:length(labels)
%     tmp_ind = strcmp(distance(:,1),labels{i});
%     tmp_label = distance(tmp_ind,2);
%     tmp_num = distance_num(tmp_ind);
%     [A B] = ismember(labels, tmp_label);
%     distMat(i,A) = tmp_num(B(A));
% end
% for i = 1:length(labels)-1
%     for j = i+1:length(labels)
%         if isnan(distMat(i,j))
%             error('!')
%         end
%     end
% end
% %reformat label
% labels = cellfun(@(x) [x(1:end-1),'_',x(end)],labels,'UniformOutput',false);
% % expand the matrix
% for i = 1:length(labels)
%     for j = 1:i-1
%         distMat(i,j) = distMat(j,i);
%     end
%     distMat(i,i) = 0;
% end
% save('distanceMatrix.mat','distMat','labels');

%% run flux efficiency
% prepare the model :: the input model should be constrianed model
%apply some constaints hacked from safak's file
model = changeRxnBounds(model,'EXC0050',-1000,'l');
model = changeRxnBounds(model,'EXC0050',-1000,'l');
model = changeRxnBounds(model,'RM04432',0,'l');
model = changeRxnBounds(model,'RM00112',1000,'u');
model = changeRxnBounds(model,'TCE5092',-1000,'l');
model = changeRxnBounds(model,'TCE5075',-1000,'l');
model = changeRxnBounds(model,'RM00112',-1000,'l');
model = changeRxnBounds(model,'EXC0152',0,'l');
bacMW=966.28583751;
model.S(end, strcmp('EXC0050',model.rxns)) = 0.02*bacMW*0.01; %the side cap
model.S(ismember(model.mets,{'atp[c]','h2o[c]','adp[c]','h[c]','pi[c]'}), strcmp('DGR0007',model.rxns)) = 0; 



%remove the NGAM 
model = changeRxnBounds(model,'RCC0005',0,'l');
%remove the GAM
model = addReaction(model,'BIO0010','reactionFormula','Biomass[c] ->','geneRule', 'NA','printLevel',0);
%decouple the specific biomass component by drianing them
for mymet = biomassList
    model = addDemandReaction(model,mymet,0);
end
maxDist = 25;
n = 1;
targetRxns = {'RM04432'};

changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);

[fluxEfficiency,fluxEfficiency_plus] = FluxEfficiency(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist);

