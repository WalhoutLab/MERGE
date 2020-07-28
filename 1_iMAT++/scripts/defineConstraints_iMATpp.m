function model = defineConstraints_iMATpp(model, infDefault,smallFluxDefault)
% This function is to define the uptake constraints for a native human
% model RECON2.2. It is not designed or tested for any other model.
%
% USAGE:
%
%    model = defineConstriants(model, infDefault,smallFluxDefault)
%
% INPUTS:
%    model:             input RECON2.2 model (COBRA model structure)
%    infDefault:        the default value for infinite fluxes
%    smallFluxDefault:  the default value for trace uptake fluxes
%
% OUTPUT:
%   model:              the constrained model
%
% `Yilmaz et al. (2020). Final Tittle and journal.
% .. Author: - Xuhang Li, Mar 2020

% since we are not aiming at precisely defining nutrient condition in human
% blood circulaton, we use a set of artificial constraints. The constriants
% are originally used to mimic DMEM culture media. 

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
model.lb(ismember(model.rxns,vitamins)) = -0.005;% artificially set as a small number
% set the maintenance
model = changeRxnBounds(model,'DM_atp_c_','l',1.833);% 1.833�mmol /gDW/h (Kilburn et�al., 1969)

AA = {'EX_his_L(e)';'EX_ala_L(e)';'EX_arg_L(e)';'EX_asn_L(e)';'EX_asp_L(e)';'EX_thr_L(e)';'EX_gln_L(e)';'EX_glu_L(e)';'EX_gly(e)';'EX_ile_L(e)';'EX_leu_L(e)';'EX_lys_L(e)';'EX_met_L(e)';'EX_phe_L(e)';'EX_pro_L(e)';'EX_ser_L(e)';'EX_trp_L(e)';'EX_tyr_L(e)';'EX_val_L(e)';'EX_cys_L(e)'};
model.lb(ismember(model.rxns,AA)) = -0.05;
model.lb(ismember(model.rxns,{'EX_gln_L(e)'})) = -0.5;
model.lb(ismember(model.rxns,{'EX_gthrd(e)'})) = -0.05;
model.lb(ismember(model.rxns,{'EX_glc(e)'})) = -5;% major carbon source in the media

end