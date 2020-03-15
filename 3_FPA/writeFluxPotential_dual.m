addpath /share/pkg/gurobi/900/linux64/matlab/
addpath ~/cobratoolbox/
addpath ./input/
addpath ./scripts/
addpath ./../input/
addpath ./../bins/
initCobraToolbox(false)
% target rxns
FPAtbl = readtable('relativeFluxPotentials.tsv','FileType','text','ReadRowNames',true);
alltarget = unique(cellfun(@(x) x(1:end-1),FPAtbl.Properties.RowNames,'UniformOutput',false));
targetExRxns = alltarget(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE|SNK|DMN)','once')),alltarget));
targetMets = alltarget(cellfun(@(x) ~isempty(regexp(x,'\[.\]','once')),alltarget));
targetMets2 = cellfun(@(x) ['NewMet_',x],targetMets,'UniformOutput',false);
%merge Ex rxns with Met (give a marker)
targetExRxns = [targetExRxns;targetMets2];
targetRxns = setdiff(alltarget,union(targetMets,targetExRxns));
%%
% prepare files
load('Tissue.mat');
% the default constriants for dual model
model = changeRxnBounds(model,'EXC0050_L',-1000,'l');
model = changeRxnBounds(model,'EX00001_E',-1000,'l');
model = changeRxnBounds(model,'EX00007_E',-1000,'l');
model = changeRxnBounds(model,'EX00009_E',-1000,'l');
model = changeRxnBounds(model,'EX00080_E',-1000,'l');
model = changeRxnBounds(model,'RM00112_I',0,'l');
model = changeRxnBounds(model,'RM00112_X',0,'l');
model = changeRxnBounds(model,'RM00112_I',0,'u');
model = changeRxnBounds(model,'RM00112_X',0,'u');
model = changeRxnBounds(model,'DMN0033_I',0,'l');
model = changeRxnBounds(model,'DMN0033_X',0,'l');
model = changeRxnBounds(model,'DMN0033_I',1,'u');
model = changeRxnBounds(model,'DMN0033_X',1,'u');
model = changeRxnBounds(model,'RMC0005_I',0,'l');
model = changeRxnBounds(model,'RMC0005_X',0,'l');
model = changeRxnBounds(model,'RMC0005_I',0,'u');
model = changeRxnBounds(model,'RMC0005_X',0,'u');
% parseGPR takes hugh amount of time, so preparse and integrate with model
% here 
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

% this loads the raw directional distance matrix, first need to take the
% minimum for future use
% to load from output of distance calculator, uncomment the following line
% distance_raw = readtable('distance_fromRxnToRxn.tsv','FileType','text','ReadRowNames',true);
load('./input/distance_raw.mat')
labels = distance_raw.Properties.VariableNames;
labels = cellfun(@(x) [x(1:end-1),'_',x(end)],labels,'UniformOutput',false);
distMat_raw = table2array(distance_raw);
distMat_min = zeros(size(distance_raw,1),size(distance_raw,2));
for i = 1:size(distMat_min,1)
    for j = 1:size(distMat_min,2)
        distMat_min(i,j) = min([distMat_raw(i,j),distMat_raw(j,i)]);
    end
end
distMat = distMat_min;

% set the inter-tissue distance to zero
X_ind = contains(labels,'_X_');
I_ind = contains(labels,'_I_');
distMat(X_ind,I_ind) = 0;
distMat(I_ind,X_ind) = 0;

% load the special penalties
manualPenalty = table2cell(readtable('manualPenalty.csv','ReadVariableNames',false,'Delimiter',','));
% load the special distance
manualDist = table2cell(readtable('manualDist.csv','ReadVariableNames',false,'Delimiter',','));

% load the confidence level
confidenceTbl = readtable('confidenceTable.tsv','FileType','text','ReadRowNames',true);
confidenceTbl.Properties.RowNames = cellfun(@(x) [x(1:end-1),'_',x(end)],confidenceTbl.Properties.RowNames,'UniformOutput',false);

%% prepare expression files
expressionTbl = readtable('expressionTable.tsv','FileType','text','ReadRowNames',true);
tissueLabel = expressionTbl.Properties.VariableNames;
master_expression = {};
for i = 1:length(tissueLabel)
    expression = struct();
    % NOTE: the expression (penalty) of _I compartment rxns is handled
    % specially by the constantPenalty
    expression.genes = cellfun(@(x) [x,'_X'], expressionTbl.Properties.RowNames, 'UniformOutput', false);
    expression.value = expressionTbl.(tissueLabel{i});
    master_expression{i} = expression;
end
%% run flux efficiency
% prepare the model :: the input model should be constrianed model
model = changeRxnBounds(model,'RM00112_X',1000,'u');
model = changeRxnBounds(model,'RM00112_I',1000,'u');
model = changeRxnBounds(model,'RM00112_X',-1000,'l');
model = changeRxnBounds(model,'RM00112_I',-1000,'l');

bacMW=966.28583751;
model.S(end, strcmp('EXC0050_L',model.rxns)) = 0.02*bacMW*0.01; %the side cap
model.S(ismember(model.mets,{'atp_I[c]','h2o_I[c]','adp_I[c]','h_I[c]','pi_I[c]'}), strcmp('DGR0007_L',model.rxns)) = 0; 

%remove the NGAM 
model = changeRxnBounds(model,'RCC0005_X',0,'l');
model = changeRxnBounds(model,'RCC0005_I',0,'l');
%% block the -2/-3 reactions
model_irrev = convertToIrreversible(model);
%unify naminclature
tmp_ind = ~cellfun(@(x) any(regexp(x,'_(f|b|r)$')),model_irrev.rxns); %has no suffix
model_irrev.rxns(tmp_ind) = cellfun(@(x) [x,'_f'],model_irrev.rxns(tmp_ind), 'UniformOutput',false);
model_irrev.rxns = regexprep(model_irrev.rxns, '_b$','_r');
X_rxns = model_irrev.rxns(contains(model_irrev.rxns,'_X_'));
for i = 1:length(tissueLabel)
    blockList{i} = confidenceTbl.Properties.RowNames((confidenceTbl.(tissueLabel{i}) <= -2));
    if strcmp(tissueLabel{i},'Intestine')
        blockList{i} = [blockList{i};X_rxns]; %block all X tissue reactions for intestine
    end     
end
blockList{end+1} = {};%for supercond, dont block anything when optimizing X tissue
blockList_I = blockList;
blockList_I{end} = X_rxns;%still block all X tissue when optimzing supercond of a Intestine tissue
%% prepare paremeters
n = 1.5;
maxDist = 31;

changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);

myCluster = parcluster('local');
myCluster.NumWorkers = 128;
saveProfile(myCluster);
parpool(39,'SpmdEnabled',false);
%% first lets do the normal rxn FPA
%% OPTIMIZE NON-INTESTINAL side
fprintf('Calculating the FPA for non-intestinal tissues...\n');
targetRxns_X = cellfun(@(x) [x,'_X'],targetRxns,'UniformOutput',false);
% for non-intestinal rxn, two special treatments are done:
%(1) the intestinal transpoters are ignored for flux sum by setting
%additional manual penalty (apply to all including supercond)
% load the manually defined penalty for Intestinal transporters (as small
% penalty for them when integrating other tissues)
manualPenalty2 = table2cell(readtable('nonIntestineManualPenalty.csv','ReadVariableNames',false,'Delimiter',','));
manualPenalty2 = [manualPenalty2;manualPenalty];
%(2) set a uniform penalty for every tissue's intestinal rxns by the
%constant penalty
% doing special treatment for intestine tissue - this is a constant penalty
% (intestinal rxn has the same penalty in all tissue optimizations; it
% could also be the reaction that we would like to keep a constant penalty
% (so they are independent of relative reaction) across conditions
% across all the optimizations 
% here we pre-calculate the intestine penalty and then overaly it in the
% calculation
% we cheat the program by a mimic expression profile to calculate intestine
% expression 
master_expression2 = {};
for i = 1:length(tissueLabel)
    expression = struct();
    % NOTE: the expression (penalty) of _I compartment rxns is handled
    % specially by the constantPenalty
    expression.genes = cellfun(@(x) [x,'_I'], expressionTbl.Properties.RowNames, 'UniformOutput', false);
    expression.value = expressionTbl.(tissueLabel{i});
    master_expression2{i} = expression;
end
penalty = calculatePenalty(model,master_expression2,manualPenalty2);
I_ind = contains(model.rxns,'_I');
constantPenalty2 = [model.rxns(I_ind),mat2cell(penalty(I_ind,strcmp(tissueLabel,'Intestine')),ones(sum(I_ind),1))];

[fluxEfficiency_X,fluxEfficiency_plus_X] = FPA(model,targetRxns_X,master_expression,distMat,labels,n, manualPenalty2,manualDist,maxDist,blockList,constantPenalty2);
%% OPTIMIZE INTESTINAL side
fprintf('Calculating the FPA for intestine...\n');
% for intestinal rxns, no special treatment is done. but to apply penalty
% (gene expression), we use the constant penalty here (the normal ones are
% used for X compartment)
% setup intestinal rxn penalty by constant penalty (cheat the program by
% constant penalty since it will overide the default penalties
targetRxns_I = cellfun(@(x) [x,'_I'],targetRxns,'UniformOutput',false);
penalty = calculatePenalty(model,master_expression2,manualPenalty);
I_ind = contains(model.rxns,'_I');
constantPenalty = [model.rxns(I_ind),mat2cell(penalty(I_ind,strcmp(tissueLabel,'Intestine')),ones(sum(I_ind),1))];

[fluxEfficiency_I,fluxEfficiency_plus_I] = FPA(model,targetRxns_I,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList_I,constantPenalty);
%% first lets do the metabolite centric FPA
% the difference is specific rxn(s) need to be blocked to prevent loop and
% false FPA value. these specific rxns are related to the objective in the
% optimization (producing or degradating); hence, we need to optimize the
% forward and/or reverse direction of a transporter seperately.
% we are mainly interested in the transportation of organic molecules. some
% co-transporting metabolite are not the target of interest so that they
% shouldnt be blocked(limiting)
helperMet = {'co2','amp','nadp','nadph','ppi','o2','nadh','nad','pi','adp','coa', 'atp', 'h2o', 'h', 'gtp', 'gdp','etfrd','etfox','na1','oh1','cl','i'}; %define a set of co-transported metabolite that are not the metabolite of interest
penalty_defined_X = calculatePenalty(model,master_expression,manualPenalty2);
penalty_defined_I = calculatePenalty(model,master_expression,manualPenalty);
fluxEfficiency_X_met = cell(length(targetExRxns),size(penalty,2));
fluxEfficiency_plus_X_met = cell(length(targetExRxns),size(penalty,2));
fluxEfficiency_I_met = cell(length(targetExRxns),size(penalty,2));
fluxEfficiency_plus_I_met = cell(length(targetExRxns),size(penalty,2));
environment = getEnvironment();
parfor i = 1:length(targetExRxns)
    restoreEnvironment(environment);
    %% demand reaction
    if contains(targetExRxns{i},'DMN') || contains(targetExRxns{i},'SNK') %demand reaction; note in the dual model, sink is actually demand, the uptake side is replaced by UPK
        % optimize the forward (producing potential)
        % all the producing end reaction
        % find the metabolite being transported/drained - the reactant for
        % demand reaction
        myMet = setdiff(regexprep(model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_X'])) < 0),'_.\[.\]',''),helperMet);
        % find the associated end reactions (producing the metabolite)
        myRxns = {};
        for j = 1:length(myMet)
            myRxns = [myRxns;model.rxns(any(model.S(cellfun(@(x) ~isempty(regexp(x,['^',myMet{j},'_.\[.\]$'],'once')),model.mets),:),1))]; %all rxns use these mets
        end
        myRxns = myRxns(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE)','once')),myRxns));
        % block these reactions and optimize 
        % only block the reactions for corresponding tissue
        model_tmp = model;
        myRxns_X = myRxns(cellfun(@(x) ~isempty(regexp(x,'(_X)$','once')),myRxns));
        model_tmp.ub(ismember(model_tmp.rxns,myRxns_X)) = 0;
        model_tmp.lb(ismember(model_tmp.rxns,myRxns_X)) = 0;
        [fluxEfficiency_X_met(i,:),fluxEfficiency_plus_X_met(i,:)] = FPA(model_tmp,{[targetExRxns{i},'_X']},master_expression,distMat,labels,n, manualPenalty2,manualDist,maxDist,blockList,constantPenalty2,false,penalty_defined_X);
        model_tmp = model;
        myRxns_I = myRxns(cellfun(@(x) ~isempty(regexp(x,'(_L|_I)$','once')),myRxns));
        model_tmp.ub(ismember(model_tmp.rxns,myRxns_I)) = 0;
        model_tmp.lb(ismember(model_tmp.rxns,myRxns_I)) = 0;
        [fluxEfficiency_I_met(i,:),fluxEfficiency_plus_I_met(i,:)] = FPA(model_tmp,{[targetExRxns{i},'_I']},master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList_I,constantPenalty,false,penalty_defined_I);
    %% sink or uptake reaction
    elseif contains(targetExRxns{i},'UPK') %sink or uptake reaction
        % note we set the distance for all uptake rxn to be zero. we need
        % to remove this manual distance for a uptake rxn when this rxn is
        % in question
        manualDist_tmp_X = manualDist(~strcmp(manualDist(:,1),[targetExRxns{i},'_X_r']),:);
        manualDist_tmp_I = manualDist(~strcmp(manualDist(:,1),[targetExRxns{i},'_I_r']),:);
        % optimize the reverse (degradation potential)
        % all the degradation end reaction (secretion)
        % find the metabolite being transported/drained - also reactant 
        myMet = setdiff(regexprep(model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_X']))<0),'_.\[.\]',''),helperMet);
        % the pseudo-mets like "sideMet" are not removed, since their
        % associated reactions are uptakes that wouldnt be constrained 
        % find the associated end reactions (producing the metabolite)
        myRxns = {};
        for j = 1:length(myMet)
            metInd = cellfun(@(x) ~isempty(regexp(x,['^',myMet{j},'_.\[.\]$'],'once')),model.mets);
            rxnCandidate = model.rxns(any(model.S(metInd,:),1));
            %some very special SNK and DMN contains the metabolite on the
            %right side, which is not actually draining the metabolite. ==>
            %kick them out
            rxnCandidate(cellfun(@(x) ~isempty(regexp(x,'^(DMN|SNK)','once')),rxnCandidate) & ~any(model.S(metInd,ismember(model.rxns,rxnCandidate)) < 0,1)') = [];
            myRxns = [myRxns;rxnCandidate]; %all rxns use these mets           
        end
        myRxns = myRxns(cellfun(@(x) ~isempty(regexp(x,'^(DMN|TCE|SNK)','once')),myRxns));
        % sinkPatch = myRxns(cellfun(@(x) ~isempty(regexp(x,'^SNK','once')),myRxns));
        % block these reactions and optimize 
        model_tmp = model;
        myRxns_X = myRxns(cellfun(@(x) ~isempty(regexp(x,'(_X)$','once')),myRxns));
        model_tmp.ub(ismember(model_tmp.rxns,myRxns_X)) = 0;
        model_tmp.lb(ismember(model_tmp.rxns,myRxns_X)) = 0;
        [fluxEfficiency_X_met(i,:),fluxEfficiency_plus_X_met(i,:)] = FPA(model_tmp,{[targetExRxns{i},'_X']},master_expression,distMat,labels,n, manualPenalty2,manualDist_tmp_X,maxDist,blockList,constantPenalty2,false,penalty_defined_X);
        model_tmp = model;
        myRxns_I = myRxns(cellfun(@(x) ~isempty(regexp(x,'(_L|_I)$','once')),myRxns));
        model_tmp.ub(ismember(model_tmp.rxns,myRxns_I)) = 0;
        model_tmp.lb(ismember(model_tmp.rxns,myRxns_I)) = 0;
        [fluxEfficiency_I_met(i,:),fluxEfficiency_plus_I_met(i,:)] = FPA(model_tmp,{[targetExRxns{i},'_I']},master_expression,distMat,labels,n, manualPenalty,manualDist_tmp_I,maxDist,blockList_I,constantPenalty,false,penalty_defined_I);
    %% transporter reactions
    elseif contains(targetExRxns{i},'TCE') %transporter reactions
        % transporters are complicated as they have two directions and two
        % tissue locations; we calculate these four conditions individually
        % the second problem is that multiple metabolites might be
        % transported in different directions, so we need to treat them
        % seprately
        % find the metabolite being transported/drained
        myMetSet = setdiff(regexprep(model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_X']))<0),'_.\[.\]',''),helperMet);
        %% forward direction potential
        % first obtain the X tissue potential  
        % find the associated end reactions (producing the metabolite)
        % analyze the forward direction is exporting or importing for EACH
        % mets involved, and block the corresponding rxns
        model_tmp = model;
        targetrxn_fullName = [targetExRxns{i},'_X'];
        for k = 1:length(myMetSet) %will be skipped if no myMet was found
            myMet = myMetSet{k};
            %determine if the metabolite is exported or imported
            isExport = ismember({[myMet,'_X[c]']},model.mets(model.S(:,strcmp(model.rxns,targetrxn_fullName))<0)); %if the left side of the reaction is the [c] labeled met
            %find all associated rxns
            metInd = cellfun(@(x) ~isempty(regexp(x,['^',myMet,'_.\[.\]$'],'once')),model.mets);
            myRxns = model.rxns(any(model.S(metInd,:),1));
            %some very special SNK and DMN contains the metabolite on the
            %right side, which is not actually draining the metabolite. ==>
            %kick them out
            myRxns(cellfun(@(x) ~isempty(regexp(x,'^(DMN|SNK)','once')),myRxns) & ~any(model.S(metInd,ismember(model.rxns,myRxns)) < 0,1)') = [];
            if isExport %metabolite is exported
                ToBlock = setdiff(myRxns(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE)','once')),myRxns)),...
                    {targetrxn_fullName});%block import (except for the target TCE itself)
                % only block for X tissue
                ToBlock = ToBlock(cellfun(@(x) ~isempty(regexp(x,'(_X)$','once')),ToBlock));
                model_tmp.ub(ismember(model_tmp.rxns,ToBlock)) = 0;
                model_tmp.lb(ismember(model_tmp.rxns,ToBlock)) = 0;
            else
                ToBlock = setdiff(myRxns(cellfun(@(x) ~isempty(regexp(x,'^(DMN|TCE|SNK)','once')),myRxns)),...
                    {targetrxn_fullName});%block export (except for the target TCE itself)
                %sinkPatch = myRxns(cellfun(@(x) ~isempty(regexp(x,'^SNK','once')),myRxns));
                ToBlock = ToBlock(cellfun(@(x) ~isempty(regexp(x,'(_X)$','once')),ToBlock));
                %sinkPatch = sinkPatch(cellfun(@(x) ~isempty(regexp(x,'(_X)$','once')),sinkPatch));
                model_tmp.ub(ismember(model_tmp.rxns,ToBlock)) = 0;
                model_tmp.lb(ismember(model_tmp.rxns,ToBlock)) = 0;
                %model_tmp.ub(ismember(model_tmp.rxns,sinkPatch)) = 0;
            end
        end
        [fluxEfficiency_X_met_f,fluxEfficiency_plus_X_met_f] = FPA(model_tmp,{targetrxn_fullName},master_expression,distMat,labels,n, manualPenalty2,manualDist,maxDist,blockList,constantPenalty2,false,penalty_defined_X);

        % then obtain the I tissue potential  
        % first to determine which reaction to target, lumen side or
        % extracellular side
        if any(strcmp([targetExRxns{i},'_L'],model.rxns)) % lumen transporter exists
            if length(myMetSet)==1 %only one metabolite is tranported, so we know the target met
                %check if it is being imported
                if ismember({[myMetSet{1},'_I[c]']},model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_L']))>0)) %if the right side of the reaction is the [c] labeled met
                    targetrxn_fullName = [targetExRxns{i},'_L']; %use the lumen side for the forward direction
                else
                    targetrxn_fullName = [targetExRxns{i},'_I']; %use the EX side for the forward direction (secreting)
                end
            elseif isempty(myMetSet) %one or more minor met is involved. now see if the target met can be determined
                if length(regexprep(model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_X']))<0),'_.\[.\]','')) == 1 %only one involved metabolite
                    tmpMet = regexprep(model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_X']))<0),'_.\[.\]','');
                    %check if it is being imported
                    if ismember({[tmpMet{:},'_I[c]']},model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_L']))>0)) %if the right side of the reaction is the [c] labeled met
                        targetrxn_fullName = [targetExRxns{i},'_L']; %use the lumen side for the forward direction
                    else
                        targetrxn_fullName = [targetExRxns{i},'_I']; %use the EX side for the forward direction (secreting)
                    end
                else
                    targetrxn_fullName = [targetExRxns{i},'_I'];%cannot decide, use EX side
                end
            else
                targetrxn_fullName = [targetExRxns{i},'_I']; %use the EX side for undetermined met target
            end
        else
            targetrxn_fullName = [targetExRxns{i},'_I'];%still calculate it even though it should be zero
        end
        model_tmp = model;
        for k = 1:length(myMetSet) %will be skipped if no myMet was found
            myMet = myMetSet{k};
            %determine if the metabolite is exported or imported
            isExport = ismember({[myMet,'_I[c]']},model.mets(model.S(:,strcmp(model.rxns,targetrxn_fullName))<0)); %if the left side of the reaction is the [e] labeled met
            %find all associated rxns
            metInd = cellfun(@(x) ~isempty(regexp(x,['^',myMet,'_.\[.\]$'],'once')),model.mets);
            myRxns = model.rxns(any(model.S(metInd,:),1));
            %some very special SNK and DMN contains the metabolite on the
            %right side, which is not actually draining the metabolite. ==>
            %kick them out
            myRxns(cellfun(@(x) ~isempty(regexp(x,'^(DMN|SNK)','once')),myRxns) & ~any(model.S(metInd,ismember(model.rxns,myRxns)) < 0,1)') = [];
            if isExport %metabolite is exported
                ToBlock = setdiff(myRxns(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE)','once')),myRxns)),...
                    {targetrxn_fullName});%block import (except for the target TCE itself)
                % only block for I tissue
                ToBlock = ToBlock(cellfun(@(x) ~isempty(regexp(x,'(_I|_L)$','once')),ToBlock));
                model_tmp.ub(ismember(model_tmp.rxns,ToBlock)) = 0;
                model_tmp.lb(ismember(model_tmp.rxns,ToBlock)) = 0;
            else
                ToBlock = setdiff(myRxns(cellfun(@(x) ~isempty(regexp(x,'^(DMN|TCE|SNK)','once')),myRxns)),...
                    {targetrxn_fullName});%block export (except for the target TCE itself)
                % sinkPatch = myRxns(cellfun(@(x) ~isempty(regexp(x,'^SNK','once')),myRxns));
                ToBlock = ToBlock(cellfun(@(x) ~isempty(regexp(x,'(_I|_L)$','once')),ToBlock));
                % sinkPatch = sinkPatch(cellfun(@(x) ~isempty(regexp(x,'(_I|_L)$','once')),sinkPatch));
                model_tmp.ub(ismember(model_tmp.rxns,ToBlock)) = 0;
                model_tmp.lb(ismember(model_tmp.rxns,ToBlock)) = 0;
                % model_tmp.ub(ismember(model_tmp.rxns,sinkPatch)) = 0;
            end
        end
        [fluxEfficiency_I_met_f,fluxEfficiency_plus_I_met_f] = FPA(model_tmp,{targetrxn_fullName},master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList_I,constantPenalty,false,penalty_defined_I);
        %% reverse direction potential
        model_tmp = model;
        targetrxn_fullName = [targetExRxns{i},'_X'];
        for k = 1:length(myMetSet) %will be skipped if no myMet was found
            myMet = myMetSet{k};
            %determine if the metabolite is exported or imported
            isExport = ismember({[myMet,'_X[c]']},model.mets(model.S(:,strcmp(model.rxns,targetrxn_fullName))>0)); %if the RIGHT side of the reaction is the [c] labeled met
            %find all associated rxns
            metInd = cellfun(@(x) ~isempty(regexp(x,['^',myMet,'_.\[.\]$'],'once')),model.mets);
            myRxns = model.rxns(any(model.S(metInd,:),1));
            %some very special SNK and DMN contains the metabolite on the
            %right side, which is not actually draining the metabolite. ==>
            %kick them out
            myRxns(cellfun(@(x) ~isempty(regexp(x,'^(DMN|SNK)','once')),myRxns) & ~any(model.S(metInd,ismember(model.rxns,myRxns)) < 0,1)') = [];
            if isExport %metabolite is exported
                ToBlock = setdiff(myRxns(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE)','once')),myRxns)),...
                    {targetrxn_fullName});%block import (except for the target TCE itself)
                % only block for X tissue
                ToBlock = ToBlock(cellfun(@(x) ~isempty(regexp(x,'(_X)$','once')),ToBlock));
                model_tmp.ub(ismember(model_tmp.rxns,ToBlock)) = 0;
                model_tmp.lb(ismember(model_tmp.rxns,ToBlock)) = 0;
            else
                ToBlock = setdiff(myRxns(cellfun(@(x) ~isempty(regexp(x,'^(DMN|TCE|SNK)','once')),myRxns)),...
                    {targetrxn_fullName});%block export (except for the target TCE itself)
                % sinkPatch = myRxns(cellfun(@(x) ~isempty(regexp(x,'^SNK','once')),myRxns));
                ToBlock = ToBlock(cellfun(@(x) ~isempty(regexp(x,'(_X)$','once')),ToBlock));
                % sinkPatch = sinkPatch(cellfun(@(x) ~isempty(regexp(x,'(_X)$','once')),sinkPatch));
                model_tmp.ub(ismember(model_tmp.rxns,ToBlock)) = 0;
                model_tmp.lb(ismember(model_tmp.rxns,ToBlock)) = 0;
                % model_tmp.ub(ismember(model_tmp.rxns,sinkPatch)) = 0;
            end
        end
        [fluxEfficiency_X_met_r,fluxEfficiency_plus_X_met_r] = FPA(model_tmp,{targetrxn_fullName},master_expression,distMat,labels,n, manualPenalty2,manualDist,maxDist,blockList,constantPenalty2,false,penalty_defined_X);

        % then obtain the I tissue potential  
        % first to determine which reaction to target, lumen side or
        % extracellular side
        if any(strcmp([targetExRxns{i},'_L'],model.rxns)) % lumen transporter exists
            if length(myMetSet)==1 %only one metabolite is tranported, so we know the target met
                %check if it is being imported
                if ismember({[myMetSet{1},'_I[c]']},model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_L']))<0)) %if the left side of the reaction is the [c] labeled met
                    targetrxn_fullName = [targetExRxns{i},'_L']; %use the lumen side for the reverse direction
                else
                    targetrxn_fullName = [targetExRxns{i},'_I']; %use the EX side for the reverse direction (secreting)
                end
            elseif isempty(myMetSet) %one or more minor met is involved. now see if the target met can be determined
                if length(regexprep(model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_X']))<0),'_.\[.\]','')) == 1 %only one involved metabolite
                    tmpMet = regexprep(model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_X']))<0),'_.\[.\]','');
                    %check if it is being imported
                    if ismember({[tmpMet{:},'_I[c]']},model.mets(model.S(:,strcmp(model.rxns,[targetExRxns{i},'_L']))>0)) %if the right side of the reaction is the [c] labeled met
                        targetrxn_fullName = [targetExRxns{i},'_L']; %use the lumen side for the forward direction
                    else
                        targetrxn_fullName = [targetExRxns{i},'_I']; %use the EX side for the forward direction (secreting)
                    end
                end
            else
                targetrxn_fullName = [targetExRxns{i},'_I']; %use the EX side for undetermined met target
            end
        else
            targetrxn_fullName = [targetExRxns{i},'_I'];%default
        end
        model_tmp = model;
        for k = 1:length(myMetSet) %will be skipped if no myMet was found
            myMet = myMetSet{k};
            %determine if the metabolite is exported or imported
            isExport = ismember({[myMet,'_I[c]']},model.mets(model.S(:,strcmp(model.rxns,targetrxn_fullName))>0)); %if the RIGHT side of the reaction is the [e] labeled met
            %find all associated rxns
            metInd = cellfun(@(x) ~isempty(regexp(x,['^',myMet,'_.\[.\]$'],'once')),model.mets);
            myRxns = model.rxns(any(model.S(metInd,:),1));
            %some very special SNK and DMN contains the metabolite on the
            %right side, which is not actually draining the metabolite. ==>
            %kick them out
            myRxns(cellfun(@(x) ~isempty(regexp(x,'^(DMN|SNK)','once')),myRxns) & ~any(model.S(metInd,ismember(model.rxns,myRxns)) < 0,1)') = [];
            if isExport %metabolite is exported
                ToBlock = setdiff(myRxns(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE)','once')),myRxns)),...
                    {targetrxn_fullName});%block import (except for the target TCE itself)
                % only block for I tissue
                ToBlock = ToBlock(cellfun(@(x) ~isempty(regexp(x,'(_I|_L)$','once')),ToBlock));
                model_tmp.ub(ismember(model_tmp.rxns,ToBlock)) = 0;
                model_tmp.lb(ismember(model_tmp.rxns,ToBlock)) = 0;
            else
                ToBlock = setdiff(myRxns(cellfun(@(x) ~isempty(regexp(x,'^(DMN|TCE|SNK)','once')),myRxns)),...
                    {targetrxn_fullName});%block export (except for the target TCE itself)
                % sinkPatch = myRxns(cellfun(@(x) ~isempty(regexp(x,'^SNK','once')),myRxns));
                ToBlock = ToBlock(cellfun(@(x) ~isempty(regexp(x,'(_I|_L)$','once')),ToBlock));
                % sinkPatch = sinkPatch(cellfun(@(x) ~isempty(regexp(x,'(_I|_L)$','once')),sinkPatch));
                model_tmp.ub(ismember(model_tmp.rxns,ToBlock)) = 0;
                model_tmp.lb(ismember(model_tmp.rxns,ToBlock)) = 0;
                % model_tmp.ub(ismember(model_tmp.rxns,sinkPatch)) = 0;
            end
        end
        [fluxEfficiency_I_met_r,fluxEfficiency_plus_I_met_r] = FPA(model_tmp,{targetrxn_fullName},master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList_I,constantPenalty,false,penalty_defined_I);
        %% merge the potentials
        fluxEfficiency_X_tmp = cell(1,size(penalty,2));
        fluxEfficiency_plus_X_tmp = cell(1,size(penalty,2));
        fluxEfficiency_I_tmp = cell(1,size(penalty,2));
        fluxEfficiency_plus_I_tmp = cell(1,size(penalty,2));
        for z = 1:length(fluxEfficiency_X_tmp)
            fluxEfficiency_X_tmp{1,z} = [fluxEfficiency_X_met_f{z}(1),fluxEfficiency_X_met_r{z}(2)];
            fluxEfficiency_I_tmp{1,z} = [fluxEfficiency_I_met_f{z}(1),fluxEfficiency_I_met_r{z}(2)];
            fluxEfficiency_plus_X_tmp{1,z} = [fluxEfficiency_plus_X_met_f{z}(1),fluxEfficiency_plus_X_met_r{z}(2)];
            fluxEfficiency_plus_I_tmp{1,z} = [fluxEfficiency_plus_I_met_f{z}(1),fluxEfficiency_plus_I_met_r{z}(2)];
        end       
        [fluxEfficiency_X_met(i,:),fluxEfficiency_plus_X_met(i,:)] = deal(fluxEfficiency_X_tmp,fluxEfficiency_plus_X_tmp);
        [fluxEfficiency_I_met(i,:),fluxEfficiency_plus_I_met(i,:)] = deal(fluxEfficiency_I_tmp,fluxEfficiency_plus_I_tmp);
    elseif contains(targetExRxns{i},'NewMet_') %new reaction needed
        % it is like a demand reaction except for adding a new reaction
        myMet = regexprep(targetExRxns{i},'^NewMet_','');
        % find the associated end reactions (producing the metabolite)
        myRxns = model.rxns(any(model.S(cellfun(@(x) ~isempty(regexp(x,['^',myMet(1:end-3),'_.\[.\]$'],'once')),model.mets),:),1)); %all rxns use these mets
        myRxns = myRxns(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE)','once')),myRxns));
        % block these reactions and optimize 
        % only block the reactions for corresponding tissue
        model_tmp = model;
        % add the new demand 
        model_tmp = addReaction(model_tmp,['DMN',myMet,'_X'],'reactionFormula',[myMet(1:end-3),'_X',myMet(end-2:end),' ->'],'geneRule', 'NA','printLevel',0);
        targetrxn_fullName = ['DMN',myMet,'_X'];
        % the distance of this new demand rxn needs to be calculated; Note
        % that the distance is actually 1+min(distance(any rxn that
        % contains this met)
        [distMat_raw_2,labels_X_2] = updateDistance(targetrxn_fullName,model_tmp,distMat_raw,labels,maxDist);
        % take the min and block intestine
        distMat_X_2 = distMat_raw_2;
        for ii = 1:size(distMat_X_2,1)
            for jj = 1:size(distMat_X_2,2)
                distMat_X_2(ii,jj) = min([distMat_raw_2(ii,jj),distMat_raw_2(jj,ii)]);
            end
        end
        % set the inter-tissue distance to zero
        X_ind = contains(labels_X_2,'_X_');
        I_ind = contains(labels_X_2,'_I_');
        distMat_X_2(X_ind,I_ind) = 0;
        distMat_X_2(I_ind,X_ind) = 0;

        penalty_defined_X_2 = calculatePenalty(model_tmp,master_expression,manualPenalty2);%recaluclate the penalty mat because of adding reaction
        myRxns_X = myRxns(cellfun(@(x) ~isempty(regexp(x,'(_X)$','once')),myRxns));
        model_tmp.ub(ismember(model_tmp.rxns,myRxns_X)) = 0;
        model_tmp.lb(ismember(model_tmp.rxns,myRxns_X)) = 0;
        [fluxEfficiency_X_met(i,:),fluxEfficiency_plus_X_met(i,:)] = FPA(model_tmp,{targetrxn_fullName},master_expression,distMat_X_2,labels_X_2,n, manualPenalty2,manualDist,maxDist,blockList,constantPenalty2,false,penalty_defined_X_2);
        
        model_tmp = model;
        model_tmp = addReaction(model_tmp,['DMN',myMet,'_I'],'reactionFormula',[myMet(1:end-3),'_I',myMet(end-2:end),' ->'],'geneRule', 'NA','printLevel',0);
        targetrxn_fullName = ['DMN',myMet,'_I'];
        [distMat_raw_2,labels_I_2] = updateDistance(targetrxn_fullName,model_tmp,distMat_raw,labels,maxDist);
        % take the min and block intestine
        distMat_I_2 = distMat_raw_2;
        for ii = 1:size(distMat_I_2,1)
            for jj = 1:size(distMat_I_2,2)
                distMat_I_2(ii,jj) = min([distMat_raw_2(ii,jj),distMat_raw_2(jj,ii)]);
            end
        end
        % set the inter-tissue distance to zero
        X_ind = contains(labels_I_2,'_X_');
        I_ind = contains(labels_I_2,'_I_');
        distMat_I_2(X_ind,I_ind) = 0;
        distMat_I_2(I_ind,X_ind) = 0;
        penalty_defined_I_2 = calculatePenalty(model_tmp,master_expression,manualPenalty);%recaluclate the penalty
        myRxns_I = myRxns(cellfun(@(x) ~isempty(regexp(x,'(_L|_I)$','once')),myRxns));
        model_tmp.ub(ismember(model_tmp.rxns,myRxns_I)) = 0;
        model_tmp.lb(ismember(model_tmp.rxns,myRxns_I)) = 0;
        [fluxEfficiency_I_met(i,:),fluxEfficiency_plus_I_met(i,:)] = FPA(model_tmp,{targetrxn_fullName},master_expression,distMat_I_2,labels_I_2,n, manualPenalty,manualDist,maxDist,blockList_I,constantPenalty,false,penalty_defined_I_2);
    else
        error('Not a supported RXN type!');
    end
end
%% merge with the rxn-centered result
fluxEfficiency_X = [fluxEfficiency_X;fluxEfficiency_X_met];
fluxEfficiency_I = [fluxEfficiency_I;fluxEfficiency_I_met];

%% get the final result of reltaive flux potential
relFP_f = nan(size(fluxEfficiency_X,1),length(tissueLabel));%flux potential for forward rxns
relFP_r = nan(size(fluxEfficiency_X,1),length(tissueLabel));%flux potential for reverse rxns
for i = 1:size(fluxEfficiency_X,1)
    for j = 1:length(tissueLabel)
        if ~strcmp(tissueLabel{j},'Intestine')
            relFP_f(i,j) = fluxEfficiency_X{i,j}(1) ./ fluxEfficiency_X{i,end}(1);
            relFP_r(i,j) = fluxEfficiency_X{i,j}(2) ./ fluxEfficiency_X{i,end}(2);
        else
            relFP_f(i,j) = fluxEfficiency_I{i,j}(1) ./ fluxEfficiency_I{i,end}(1);
            relFP_r(i,j) = fluxEfficiency_I{i,j}(2) ./ fluxEfficiency_I{i,end}(2);
        end
    end
end
save('FPA_rxns.mat','relFP_f','relFP_r');
save('FPA_workspace.mat');
