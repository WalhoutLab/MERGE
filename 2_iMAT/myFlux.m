%% load model
addpath /opt/cobra/cobratoolbox/
addpath ~/conjoinedModel/new_matlab_functions/
addpath ~/cobratoolbox/
addpath /share/pkg/gurobi/810/linux64/matlab/
addpath ~/new_matlab_functions/
addpath /Users/xuhangli/Desktop/Walhout_Lab/conjoined_model_project/new_matlab_functions
initCobraToolbox;
load('Tissue.mat');
% worm = addDefaultConstraint(model,'default');
% worm = changeRxnBounds(worm,'RCC0005',10,'l');
% % change the side proportion
% worm.S(end-1, strcmp('EXC0050',worm.rxns)) = 0.01*9.6641833;
% worm.S(end, strcmp('EXC0050',worm.rxns)) = 0.01*9.6641833;
% % make the epsilon vector
% [epsilon_f,epsilon_r] = makeEpsilonSeq(worm,worm.rxns,0.01,0.5);

% reset some constraints to make the model ready for integration 
% release the input constraints for integration 
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
% parseGPR takes hugh amount of time, so preparse and integrate with model
% here 
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
%% laod the gene expression data
fname = './input/geneSets.json';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
str = regexprep(str,'set(','');
str = regexprep(str,')','');
GeneExpression = jsondecode(str);
% make expression set for six tissues

% hypodermis
%add compartment label and merge
GeneExpression.Hypodermis.high = strcat(GeneExpression.Hypodermis.high,'_X');
GeneExpression.Hypodermis.lenient = strcat(GeneExpression.Hypodermis.lenient,'_X');
GeneExpression.Hypodermis.low = strcat(GeneExpression.Hypodermis.low,'_X');
GeneExpression.Hypodermis.zero = strcat(GeneExpression.Hypodermis.zero,'_X');
GeneExpression.Intestine.high = strcat(GeneExpression.Intestine.high,'_I');
GeneExpression.Intestine.lenient = strcat(GeneExpression.Intestine.lenient,'_I');
GeneExpression.Intestine.low = strcat(GeneExpression.Intestine.low,'_I');
GeneExpression.Intestine.zero = strcat(GeneExpression.Intestine.zero,'_I');
%merge
ExpCatag = struct();
ExpCatag.high = [GeneExpression.Hypodermis.high];
%ExpCatag.dynamic = [GeneExpression.Hypodermis.lenient;GeneExpression.Intestine.lenient];
ExpCatag.low = [GeneExpression.Hypodermis.low;GeneExpression.Intestine.low];
ExpCatag.zero = [GeneExpression.Hypodermis.zero;GeneExpression.Intestine.zero];
%all the other genes are dynamic (no constriants on) 
ExpCatag.dynamic = setdiff(model.genes, unique([ExpCatag.high;ExpCatag.low;ExpCatag.zero]));
ExpCatag_Hypodermis = ExpCatag;

% Gonad
%add compartment label and merge
GeneExpression.Gonad.high = strcat(GeneExpression.Gonad.high,'_X');
GeneExpression.Gonad.lenient = strcat(GeneExpression.Gonad.lenient,'_X');
GeneExpression.Gonad.low = strcat(GeneExpression.Gonad.low,'_X');
GeneExpression.Gonad.zero = strcat(GeneExpression.Gonad.zero,'_X');
%merge
ExpCatag = struct();
ExpCatag.high = [GeneExpression.Gonad.high];
%ExpCatag.dynamic = [GeneExpression.Hypodermis.lenient;GeneExpression.Intestine.lenient];
ExpCatag.low = [GeneExpression.Gonad.low;GeneExpression.Intestine.low];
ExpCatag.zero = [GeneExpression.Gonad.zero;GeneExpression.Intestine.zero];
%all the other genes are dynamic (no constriants on) 
ExpCatag.dynamic = setdiff(model.genes, unique([ExpCatag.high;ExpCatag.low;ExpCatag.zero]));
ExpCatag_Gonad = ExpCatag;

% Glia
%add compartment label and merge
GeneExpression.Glia.high = strcat(GeneExpression.Glia.high,'_X');
GeneExpression.Glia.lenient = strcat(GeneExpression.Glia.lenient,'_X');
GeneExpression.Glia.low = strcat(GeneExpression.Glia.low,'_X');
GeneExpression.Glia.zero = strcat(GeneExpression.Glia.zero,'_X');
%merge
ExpCatag = struct();
ExpCatag.high = [GeneExpression.Glia.high];
%ExpCatag.dynamic = [GeneExpression.Hypodermis.lenient;GeneExpression.Intestine.lenient];
ExpCatag.low = [GeneExpression.Glia.low;GeneExpression.Intestine.low];
ExpCatag.zero = [GeneExpression.Glia.zero;GeneExpression.Intestine.zero];
%all the other genes are dynamic (no constriants on) 
ExpCatag.dynamic = setdiff(model.genes, unique([ExpCatag.high;ExpCatag.low;ExpCatag.zero]));
ExpCatag_Glia = ExpCatag;

% Pharynx
%add compartment label and merge
GeneExpression.Pharynx.high = strcat(GeneExpression.Pharynx.high,'_X');
GeneExpression.Pharynx.lenient = strcat(GeneExpression.Pharynx.lenient,'_X');
GeneExpression.Pharynx.low = strcat(GeneExpression.Pharynx.low,'_X');
GeneExpression.Pharynx.zero = strcat(GeneExpression.Pharynx.zero,'_X');
%merge
ExpCatag = struct();
ExpCatag.high = [GeneExpression.Pharynx.high];
%ExpCatag.dynamic = [GeneExpression.Hypodermis.lenient;GeneExpression.Intestine.lenient];
ExpCatag.low = [GeneExpression.Pharynx.low;GeneExpression.Intestine.low];
ExpCatag.zero = [GeneExpression.Pharynx.zero;GeneExpression.Intestine.zero];
%all the other genes are dynamic (no constriants on) 
ExpCatag.dynamic = setdiff(model.genes, unique([ExpCatag.high;ExpCatag.low;ExpCatag.zero]));
ExpCatag_Pharynx = ExpCatag;

% Body_wall_muscle
%add compartment label and merge
GeneExpression.Body_wall_muscle.high = strcat(GeneExpression.Body_wall_muscle.high,'_X');
GeneExpression.Body_wall_muscle.lenient = strcat(GeneExpression.Body_wall_muscle.lenient,'_X');
GeneExpression.Body_wall_muscle.low = strcat(GeneExpression.Body_wall_muscle.low,'_X');
GeneExpression.Body_wall_muscle.zero = strcat(GeneExpression.Body_wall_muscle.zero,'_X');
%merge
ExpCatag = struct();
ExpCatag.high = [GeneExpression.Body_wall_muscle.high];
%ExpCatag.dynamic = [GeneExpression.Hypodermis.lenient;GeneExpression.Intestine.lenient];
ExpCatag.low = [GeneExpression.Body_wall_muscle.low;GeneExpression.Intestine.low];
ExpCatag.zero = [GeneExpression.Body_wall_muscle.zero;GeneExpression.Intestine.zero];
%all the other genes are dynamic (no constriants on) 
ExpCatag.dynamic = setdiff(model.genes, unique([ExpCatag.high;ExpCatag.low;ExpCatag.zero]));
ExpCatag_Body_wall_muscle = ExpCatag;

% Neurons
%add compartment label and merge
GeneExpression.Neurons.high = strcat(GeneExpression.Neurons.high,'_X');
GeneExpression.Neurons.lenient = strcat(GeneExpression.Neurons.lenient,'_X');
GeneExpression.Neurons.low = strcat(GeneExpression.Neurons.low,'_X');
GeneExpression.Neurons.zero = strcat(GeneExpression.Neurons.zero,'_X');
%merge
ExpCatag = struct();
ExpCatag.high = [GeneExpression.Neurons.high];
%ExpCatag.dynamic = [GeneExpression.Hypodermis.lenient;GeneExpression.Intestine.lenient];
ExpCatag.low = [GeneExpression.Neurons.low;GeneExpression.Intestine.low];
ExpCatag.zero = [GeneExpression.Neurons.zero;GeneExpression.Intestine.zero];
%all the other genes are dynamic (no constriants on) 
ExpCatag.dynamic = setdiff(model.genes, unique([ExpCatag.high;ExpCatag.low;ExpCatag.zero]));
ExpCatag_Neurons = ExpCatag;

% Intestine
%merge
ExpCatag = struct();
ExpCatag.high = [GeneExpression.Intestine.high];
%ExpCatag.dynamic = [GeneExpression.Hypodermis.lenient;GeneExpression.Intestine.lenient];
ExpCatag.low = [GeneExpression.Intestine.low];
ExpCatag.zero = [GeneExpression.Intestine.zero];
%all the other genes are dynamic (no constriants on) 
ExpCatag.dynamic = setdiff(model.genes, unique([ExpCatag.high;ExpCatag.low;ExpCatag.zero]));
ExpCatag_Intestine = ExpCatag;

%merge all 
tissues = struct();
tissues.Hypodermis = ExpCatag_Hypodermis;
tissues.Gonad = ExpCatag_Gonad;
tissues.Glia = ExpCatag_Glia;
tissues.Pharynx = ExpCatag_Pharynx;
tissues.Body_wall_muscle = ExpCatag_Body_wall_muscle;
tissues.Neurons = ExpCatag_Neurons;
tissues.Intestine = ExpCatag_Intestine;

%% load epsilon
fname = './input/epsilon.json';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
str = regexprep(str,'u"','"');
epsilon = jsondecode(str);
fields = fieldnames(epsilon);
epsilon_f = 0.01*ones(length(model.rxns),1);
epsilon_r = 0.01*ones(length(model.rxns),1);
for i = 1: length(fields)
    if epsilon.(fields{i})(2) ~= 0
        epsilon_f(strcmp(model.rxns,fields{i})) = epsilon.(fields{i})(2);
    end
    if epsilon.(fields{i})(1) ~= 0
        epsilon_r(strcmp(model.rxns,fields{i})) = epsilon.(fields{i})(1);
    end
end
%% make OFD
names = fieldnames(tissues);
names = names(1:end-1);
for i = 1:length(names)
    fprintf('now starting to fit %s... \n',names{i});
    t1 = tic();
    ExpCatag = tissues.(names{i});
    doLatent = 0;
    storeProp = 0.01;
    SideProp = 0.02;
    ATPm = 10;
    doMinPFD = 1;
    latentCAP = 0.05;
    type = 'iMATplus';
    myCSM = struct(); %my context specific model
    [myCSM.OFD,myCSM.N_highFit,myCSM.N_zeroFit,myCSM.minLow,myCSM.minTotal,myCSM.OpenGene,myCSM.wasteDW,myCSM.HGenes,myCSM.RLNames,myCSM.latentRxn,myCSM.PFD,myCSM.Nfit_latent,myCSM.minTotal_OFD,myCSM.MILP] = autoIntegration_iMAT(model,doLatent,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCatag,doMinPFD,latentCAP,type);
    save([names{i},'.mat'],'myCSM');
    eval([names{i},' = myCSM;']);
    toc(t1)
    %OFD = fix(FluxDistribution_daf2 .* 1e7) ./ 1e7;
    % what's the atpm for each tissue?
    % why not every reaction has epsilon?
    % release some constraints?
    % what is lenient? overlap with low? 
    % need a dynamic category, otherwise these moderate expressed gene will be
    % overide by default value "-1"
    % how is undetectable genes processed? in iMAT, it has -1 score so that all
    % the reactions associated with "and" gate will be exluded from MILP
    % integration! it seems overide some real low expressed rxns (i.e, 0 and -1
    % gives -1, but should be minimized? undetected means zero?)
end
%% now do the intestine integration
% first get the average flux distribution
weights = [0.179980982;0.086522532;0.036874975;0.138793829;0.372097477000000;0.185730205];%needs to order by tissuesname
fluxSum_PFD = zeros(length(model.rxns),1);
for i = 1:length(names)
    eval(['fluxSum_PFD = ',names{i},'.PFD .* weights(i) + fluxSum_PFD;']);
end
IntestineModel_PFD = collapseX(model,'X',fluxSum_PFD);
%update GPR
parsedGPR = GPRparser_xl(IntestineModel_PFD);% Extracting GPR data from model
IntestineModel_PFD.parsedGPR = parsedGPR;
% update the epsilon list
epsilon_new_f = 0.01 * ones(length(IntestineModel_PFD.rxns),1);
epsilon_new_r = epsilon_new_f;
for i = 1: length(IntestineModel_PFD.rxns)
    if any(ismember(IntestineModel_PFD.rxns(i), model.rxns))
        epsilon_new_f(i) = epsilon_f(ismember(model.rxns,IntestineModel_PFD.rxns(i)));
        epsilon_new_r(i) = epsilon_r(ismember(model.rxns,IntestineModel_PFD.rxns(i)));
    end
end
% integration 
fprintf('now starting to fit %s... \n','Intestine');
tic()
ExpCatag = tissues.Intestine;
doLatent = 0;
storeProp = 0.01;
SideProp = 0.02;
ATPm = 10;
doMinPFD = 1;
latentCAP = 0.05;
type = 'canonical';
myCSM = struct(); %my context specific model
%PFD fitting
[~,myCSM.N_highFit,myCSM.N_zeroFit,myCSM.minLow,myCSM.minTotal,myCSM.OpenGene,~,myCSM.HGenes,myCSM.RLNames,~,myCSM.PFD,~,~,~] = autoIntegration_iMAT(IntestineModel_PFD,doLatent,storeProp,SideProp,epsilon_new_f,epsilon_new_r, ATPm, ExpCatag,doMinPFD,latentCAP,type);
%OFD fitting
save('Intestine.mat','myCSM');
toc()
%% load the reference data
fname = './input/fittingQuality.json';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
str = regexprep(str,'u"','"');
str = regexprep(str,'set(','');
str = regexprep(str,')','');
fittingQuality = jsondecode(str);
%% load the json
fname = './input/OFD.json';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
str = regexprep(str,'u"','"');
str = regexprep(str,'set(','');
str = regexprep(str,')','');
OFD_safak = jsondecode(str);
%% load the json
fname = './input/collapsedX_exampleConstraints.json';
str = fileread(fname);
for i = 1:length(str)
    if str(i) == "'"
        str(i) = '"';
    end
end
str = regexprep(str,'u"','"');
str = regexprep(str,'set(','');
str = regexprep(str,')','');
costraintsExp = jsondecode(str);