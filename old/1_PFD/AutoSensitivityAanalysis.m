function [sensitivity] = AutoSensitivityAanalysis(worm,ExpCatag,storageCoef,sideCoef,ATPm,epsilon_f,epsilon_r,capacity_f,capacity_r,clusterFlag,doFullSensitivity)
%% inputs
% model ... the model; model needs to be constrained
% storeProp,SideProp ... proportion of sides and storage molecules
% epsilon_f/r; capacity_f/r ... calculated epsilon/flux capacity label vector for every reaction
% ExpCatag ... the expressiom category of each gene; note that the gene
% name should in iCEL format
% clusterFlag: use cluster or not
addpath ~/cobratoolbox/
addpath /share/pkg/gurobi/810/linux64/matlab/
addpath ~/new_matlab_functions/
initCobraToolbox;
%% calculate the original integration
[~,N_highFit_ori,N_zeroFit_ori,minLow_ori,minTotal_ori,~,~,~,~,~,~,FluxDistribution_ori] = autoIntegration_latent(worm,false,storageCoef,sideCoef,epsilon_f,epsilon_r,[],[], ATPm, ExpCatag);
%% perturb every flux-capable rxn
tol = 1e-7; % the default flux numerical tolerance is alway 1e-7
if clusterFlag
    % load cluster environment
    loadCluster(60,'1024','large','1:30')
else
    parpool(2)
    pctRunOnAll initCobraToolbox;
end
%make the output list
sensitivity = struct('Name',{{}},'Score',[],'rxnID',{{}});
for i = 1:length(worm.rxns)
    if capacity_f(i)
        sensitivity.Name = [sensitivity.Name;{[worm.rxns{i},'_f']}];
        sensitivity.Score = [sensitivity.Score;nan];
        sensitivity.rxnID = [sensitivity.rxnID;worm.rxns(i)];
    end
    if capacity_r(i)
        sensitivity.Name = [sensitivity.Name;{[worm.rxns{i},'_r']}];
        sensitivity.Score = [sensitivity.Score;nan];
        sensitivity.rxnID = [sensitivity.rxnID;worm.rxns(i)];
    end
end
sensitivity_Name = sensitivity.Name;
sensitivity_rxnID = sensitivity.rxnID;
sensitivity_Score = sensitivity.Score;
parfor i = 1:length(sensitivity_Name)
    model_ptb = worm;
    rxnInd = find(strcmp(model_ptb.rxns,sensitivity_rxnID{i}));
    if abs(FluxDistribution_ori(rxnInd)) < tol %zero -> non-zero; give a negative score
        if doFullSensitivity %will skip negative score caculation if only to calculate OFD
            if strcmp(sensitivity_Name{i}(end-1:end),'_f')
                % if this reaction is a forward reaction
                model_ptb.lb(rxnInd) = epsilon_f(rxnInd);
                try %in case infeasible
                    [~,N_highFit_p,N_zeroFit_p,minLow_p,minTotal_p] = ...
                        autoIntegration_latent(model_ptb,false,storageCoef,sideCoef,epsilon_f,epsilon_r,[],[], ATPm, ExpCatag,doFullSensitivity);%whether do minTotal depends on whther do full sensitivity analysis
                    % scoring 
                    if N_highFit_p + N_zeroFit_p < N_highFit_ori + N_zeroFit_ori
                        sensitivity_Score(i) = -3;
                    elseif minLow_p > minLow_ori + 1e-7 %the numerical tolerance
                        sensitivity_Score(i) = -2;
                    elseif doFullSensitivity % calculate minTotal and compare minTotal
                        if minTotal_p > minTotal_ori + 1e-7
                            sensitivity_Score(i) = -1;
                        end
                    else
                        sensitivity_Score(i) = 0;
                    end
                catch
                    sensitivity_Score(i) = -3; %infeasible is also a strongest sensitivity
                    InfLog = fopen('./InfCases.txt','a+');%write the infeasible cases to a log file
                    fprintf(InfLog,'The rxn %s is infeasible after perturbation; the original flux is %.2f\n',sensitivity_Name{i},FluxDistribution_ori(rxnInd));
                    fclose(InfLog);
                end
            else %its a reverse reaction
                model_ptb.ub(rxnInd) = -epsilon_r(rxnInd);
                try
                    [~,N_highFit_p,N_zeroFit_p,minLow_p,minTotal_p] = ...
                        autoIntegration_latent(model_ptb,false,storageCoef,sideCoef,epsilon_f,epsilon_r,[],[], ATPm, ExpCatag,doFullSensitivity);
                    if N_highFit_p + N_zeroFit_p < N_highFit_ori + N_zeroFit_ori
                        sensitivity_Score(i) = -3;
                    elseif minLow_p > minLow_ori + 1e-7 %the numerical tolerance
                        sensitivity_Score(i) = -2;
                    elseif doFullSensitivity
                        if minTotal_p > minTotal_ori + 1e-7
                            sensitivity_Score(i) = -1;
                        end
                    else
                        sensitivity_Score(i) = 0;
                    end        
                catch
                    sensitivity_Score(i) = -3;
                    InfLog = fopen('./InfCases.txt','a+');%write the infeasible cases to a log file
                    fprintf(InfLog,'The rxn %s is infeasible after perturbation; the original flux is %.2f\n',sensitivity_Name{i},FluxDistribution_ori(rxnInd));
                    fclose(InfLog);
                end
            end
        else
            sensitivity_Score(i) = 0;
        end
    else %it is orginally non-zero flux
        if (strcmp(sensitivity_Name{i}(end-1:end),'_f') && FluxDistribution_ori(rxnInd) > 0) || (strcmp(sensitivity_Name{i}(end-1:end),'_r') && FluxDistribution_ori(rxnInd) < 0) %f reaction and ori flux is f || r reaction and ori flux is r
            model_ptb.lb(rxnInd) = 0;
            model_ptb.ub(rxnInd) = 0; %just put it zero
            try
                [~,N_highFit_p,N_zeroFit_p,minLow_p,minTotal_p] = ...
                    autoIntegration_latent(model_ptb,false,storageCoef,sideCoef,epsilon_f,epsilon_r,[],[], ATPm, ExpCatag,doFullSensitivity);
                if N_highFit_p + N_zeroFit_p < N_highFit_ori + N_zeroFit_ori
                    sensitivity_Score(i) = 3;
                elseif minLow_p > minLow_ori + 1e-7 %the numerical tolerance
                    sensitivity_Score(i) = 2;
                elseif doFullSensitivity
                    if minTotal_p > minTotal_ori + 1e-7
                        sensitivity_Score(i) = 1;
                    end
                else
                    sensitivity_Score(i) = 0;
                end
            catch
                sensitivity_Score(i) = 3; 
                InfLog = fopen('./InfCases.txt','a+');%write the infeasible cases to a log file
                fprintf(InfLog,'The rxn %s is infeasible after perturbation; the original flux is %.2f\n',sensitivity_Name{i},FluxDistribution_ori(rxnInd));
                fclose(InfLog);
            end
        else %is the nonflux reverse reaction, put it to the corresponding epsilon flux
            if doFullSensitivity
                if strcmp(sensitivity_Name{i}(end-1:end),'_f')
                    model_ptb.lb(rxnInd) = epsilon_f(rxnInd);
                else
                    model_ptb.ub(rxnInd) = -epsilon_r(rxnInd);
                end
                try %in case infeasible
                    [~,N_highFit_p,N_zeroFit_p,minLow_p,minTotal_p] = ...
                        autoIntegration_latent(model_ptb,false,storageCoef,sideCoef,epsilon_f,epsilon_r,[],[], ATPm, ExpCatag,doFullSensitivity);
                    % scoring 
                    if N_highFit_p + N_zeroFit_p < N_highFit_ori + N_zeroFit_ori
                        sensitivity_Score(i) = -3;
                    elseif minLow_p > minLow_ori + 1e-7 %the numerical tolerance
                        sensitivity_Score(i) = -2;
                    elseif doFullSensitivity
                        if minTotal_p > minTotal_ori + 1e-7
                            sensitivity_Score(i) = -1;
                        end
                    else
                        sensitivity_Score(i) = 0;
                    end
                catch
                    sensitivity_Score(i) = -3; %infeasible is also a strongest sensitivity
                    InfLog = fopen('./InfCases.txt','a+');%write the infeasible cases to a log file
                    fprintf(InfLog,'The rxn %s is infeasible after perturbation; the original flux is %.2f\n',sensitivity_Name{i},FluxDistribution_ori(rxnInd));
                    fclose(InfLog);
                end
            else
                sensitivity_Score(i) = 0;
            end
        end
    end        
end
sensitivity.Name = sensitivity_Name;
sensitivity.rxnID = sensitivity_rxnID;
sensitivity.Score = sensitivity_Score;
end
