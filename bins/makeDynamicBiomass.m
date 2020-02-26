%it is based on iJOBiomass but should be adapted to iCEL biomass

%in the input model, inquired biomass should be modified in a way that
%growth maintance is deleted and the right equation is only by-product of
%biomass formation
function new_model  = makeDynamicBiomass(model,target_Biomass, GM_coef)
new_model = model;
new_model = addMultipleMetabolites(new_model,{'biomass_mass','biomass_GM'});
%add the pseudoreactions that drain seperate biomass component 
LeftMetInd = find(new_model.S(:,strcmp(new_model.rxns,target_Biomass)) < 0 );
RightMetInd = find(new_model.S(:,strcmp(new_model.rxns,target_Biomass)) > 0 );
for i = LeftMetInd'
    myMW = MolMass(new_model.metFormulas{i});
    new_model = addReaction(new_model, ['biomass_L_part_',new_model.mets{i}], 'metaboliteList',{new_model.mets{i},'biomass_mass'},'stoichCoeffList',[-1 myMW/1000], 'reversible',false);
end
for i = RightMetInd'
    myMW = MolMass(new_model.metFormulas{i});
    new_model = addReaction(new_model, ['biomass_R_part_',new_model.mets{i}], 'metaboliteList',{new_model.mets{i},'biomass_mass'},'stoichCoeffList',[1 -myMW/1000], 'reversible',false);
end
%add reactions for GM part
new_model = addReaction(new_model, 'biomass_GM_comsumption', 'metaboliteList',{'atp_c','h2o_c','adp_c','pi_c','h_c','biomass_GM'},'stoichCoeffList',[-1 -1 1 1 1 1], 'reversible',false);
new_model = addReaction(new_model, 'biomass_drain', 'metaboliteList',[{'biomass_mass','biomass_GM'},new_model.mets(RightMetInd)'],'stoichCoeffList',[-1 -GM_coef repmat(1,1,length(RightMetInd))], 'reversible',false);
%block the old biomass reaction
new_model.rxns(strcmp(new_model.rxns,target_Biomass)) = {[target_Biomass,'_old']};
new_model = changeRxnBounds(new_model,[target_Biomass,'_old'],0,'l');
new_model = changeRxnBounds(new_model,[target_Biomass,'_old'],0,'u');
new_model.csense = [new_model.csense,'EE'];
end