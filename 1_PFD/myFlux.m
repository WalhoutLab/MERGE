%% load model
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
ExpCatag_hypodermis = ExpCatag;

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

%merge all 
tissues = struct();
tissues.hypodermis = ExpCatag_hypodermis;
tissues.Gonad = ExpCatag_Gonad;
tissues.Glia = ExpCatag_Glia;
tissues.Pharynx = ExpCatag_Pharynx;
tissues.Body_wall_muscle = ExpCatag_Body_wall_muscle;
tissues.Neurons = ExpCatag_Neurons;
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
for i = 1:length(names)
    fprintf('now starting to fit %s... \n',names{i});
    ExpCatag = tissues.(names{i});
    doLatent = 1;
    storeProp = 0.01;
    SideProp = 0.02;
    capacity_f = [];
    capacity_r = [];
    ATPm = 10;
    doMinPFD = 1;
    doFullSensitivity = 0;
    myCSM = struct(); %my context specific model
    [myCSM.OFD,myCSM.N_highFit,myCSM.N_zeroFit,myCSM.minLow,myCSM.minTotal,myCSM.OpenGene,myCSM.wasteDW,myCSM.HGenes,myCSM.RLNames,myCSM.latentRxn,myCSM.PFD,myCSM.Nfit_latent,myCSM.minTotal_OFD,myCSM.MILP] = autoIntegration_latent(model,doLatent,storeProp,SideProp,epsilon_f,epsilon_r,capacity_f,capacity_r, ATPm, ExpCatag,doMinPFD,doFullSensitivity);
    save([names{i},'.mat'],'myCSM');
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
