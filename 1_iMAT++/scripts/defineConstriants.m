function model = defineConstriants(model, infDefault,smallFluxDefault, MFAdata)
% This function is to define the uptake constrainst for a native human
% model RECON2.2. It is not designed or tested for any other model.
%
% USAGE:
%
%    model = defineConstriants(model, infDefault,smallFluxDefault, FVA)
%
% INPUTS:
%    model:             input RECON2.2 model (COBRA model structure)
%    infDefault:        the default value for infinite fluxes
%    smallFluxDefault:  the default value for trace uptake fluxes
%    MFA:               the Metabolic Flux Measurement data
%
% OUTPUT:
%   model:              the constrianed model
%
% `Yilmaz et al. (2020). Final Tittle and journal.
% .. Author: - Xuhang Li, Mar 2020

% set infinite
model.ub(isinf(model.ub)) = infDefault;
model.lb(isinf(model.lb)) = -infDefault;
% open all exchange to a small flux default
model.lb(cellfun(@(x) ~isempty(regexp(x,'^EX_','once')),model.rxns)) = -smallFluxDefault;
model.lb(cellfun(@(x) ~isempty(regexp(x,'^sink_','once')),model.rxns)) = -smallFluxDefault;
% define the freely avaiable inorganic media content 
media = {'EX_ca2(e)',...
        'EX_cl(e)',...
        'EX_fe2(e)',...
        'EX_fe3(e)',...
        'EX_h(e)',...
        'EX_h2o(e)',...
        'EX_k(e)',...
        'EX_na1(e)',...
        'EX_nh4(e)',...
        'EX_so4(e)',...
        'EX_pi(e)',...
        'EX_o2(e)'};
model.lb(ismember(model.rxns,media)) = -infDefault;% media ion set to free

% define the vitamin input
vitamins = {'EX_btn(e)',...
        'EX_chol(e)',...
        'EX_pnto_R(e)',...
        'EX_fol(e)',...
        'EX_ncam(e)',...
        'EX_pydxn(e)',...
        'EX_ribflv(e)',...
        'EX_thm(e)',...
        'EX_adpcbl(e)',...
        };
model.lb(ismember(model.rxns,vitamins)) = -1;%artificially set as -1
% set the maintaince
model = changeRxnBounds(model,'DM_atp_c_','l',1.833);%1.833 mmol gDW?1 h?1 (Kilburn et al., 1969)
% the flux in the analysis is reported as fmol/cell/h; accoring to https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005698#sec011
% 1mmol/gdw ~ 10^12fmol/10^12 cells ~ 1 fmol/cell
% because only cell volume data is available yet the density is not, we
% ignore the dW/cell variation across 60 lines. so the growth rate is

% we use the minimal uptake rate for a metabolite when multiple conditions are provided (this will be used to calculate epsilon)
% and the constriant will be the real uptake rate when only one condition
% is provided
for i = 1:size(MFAdata,1)
    myrxn = MFAdata.rxnID{i};
    fluxes = MFAdata{i,6:size(MFAdata,2)};
    if any(fluxes < 0) %could be uptaken
        model.lb(strcmp(model.rxns,myrxn)) = max(fluxes(fluxes < 0));
    else %don't allow uptake
        model.lb(strcmp(model.rxns,myrxn)) = 0;
    end
end
if ~any(strcmp(model.rxns,'transport_dhap'))%not modified model
    % fix some conflicts between model reconstruction and the flux data
    % dhap can only be uptaken but cannot carry influx ==> add a transport rxn
    model = addReaction(model,['transport_dhap'],'reactionFormula','dhap[e] <==> dhap[c]','geneRule', 'NA','printLevel',1);
    % EX_hom_L(e) cannot carry flux ==> add a celluar demand for this
    % it is a co-transporting circular met
    model = addDemandReaction(model,'hom_L[c]');
    % EX_sbt-d(e) can only be uptaken but cannot carry influx ==> change transport reversibility
    model.lb(strcmp(model.rxns,'SBTle')) = -infDefault;
end
% additionally, the essential aa histidine and nonessential aa cys has no data
% we give sufficient amout of histidine to make sure it is not limiting the
% system, but we keep cys at default due to lack of data
AA = {'EX_his_L(e)'};%,'EX_cys_L(e)'};
model.lb(ismember(model.rxns,AA)) = -10;
end
