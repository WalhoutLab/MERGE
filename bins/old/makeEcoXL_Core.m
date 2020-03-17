load('iJO1366.mat');
ecoli = iJO1366;
%add sink
%[num,txt,all] = xlsread('sinks_to_fix_rxn.xlsx','met_list');
%sinksRxns = all(2:end,1);
%ecoli = addSinkReactions(ecoli, sinksRxns);
%modify biomass
ecoli = changeRxnBounds(ecoli,'BIOMASS_Ec_iJO1366_WT_53p95M',0,'l');
ecoli = changeRxnBounds(ecoli,'BIOMASS_Ec_iJO1366_WT_53p95M',0,'u');
%seperate APT-GM from biomass reaction . -- 53.95 GAM estimate
ecoli.S(strcmp(ecoli.mets,'atp_c'),strcmp(ecoli.rxns,'BIOMASS_Ec_iJO1366_core_53p95M')) = ecoli.S(strcmp(ecoli.mets,'atp_c'),strcmp(ecoli.rxns,'BIOMASS_Ec_iJO1366_core_53p95M')) + 53.95;
ecoli.S(strcmp(ecoli.mets,'h2o_c'),strcmp(ecoli.rxns,'BIOMASS_Ec_iJO1366_core_53p95M')) = ecoli.S(strcmp(ecoli.mets,'h2o_c'),strcmp(ecoli.rxns,'BIOMASS_Ec_iJO1366_core_53p95M')) + 53.95;
ecoli.S(strcmp(ecoli.mets,'adp_c'),strcmp(ecoli.rxns,'BIOMASS_Ec_iJO1366_core_53p95M')) = 0;
ecoli.S(strcmp(ecoli.mets,'h_c'),strcmp(ecoli.rxns,'BIOMASS_Ec_iJO1366_core_53p95M')) = 0;
ecoli.S(strcmp(ecoli.mets,'pi_c'),strcmp(ecoli.rxns,'BIOMASS_Ec_iJO1366_core_53p95M')) = ecoli.S(strcmp(ecoli.mets,'pi_c'),strcmp(ecoli.rxns,'BIOMASS_Ec_iJO1366_core_53p95M'))-53.95;
%make dynamic biomass
ecoli = makeDynamicBiomass(ecoli,'BIOMASS_Ec_iJO1366_core_53p95M',53.95);
%%
ecoli = changeObjective(ecoli, 'biomass_drain');
FBA_dynamic_biomass(ecoli, 'biomass_drain','BIOMASS_Ec_iJO1366_core_53p95M_old',0.9,1.1,'max')