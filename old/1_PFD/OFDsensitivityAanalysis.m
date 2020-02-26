function [sensitivity] = OFDsensitivityAanalysis(worm,ExpCatag,storageCoef,sideCoef,ATPm,epsilon_f,epsilon_r,capacity_f,capacity_r,clusterFlag)
%% inputs
% model ... the model; model needs to be constrained
% storeProp,SideProp ... proportion of sides and storage molecules
% epsilon_f/r; capacity_f/r ... calculated epsilon/flux capacity label vector for every reaction
% ExpCatag ... the expressiom category of each gene; note that the gene
% name should in iCEL format
% clusterFlag: use cluster or not

% the ss Score is reassigned as:
% 3 --> high/zero expression violation
% 2 --> low expression violation
% 1.5 --> latent violation
% 1 --> minFlux violation
% 0 --> no violation
addpath ~/cobratoolbox/
addpath /share/pkg/gurobi/810/linux64/matlab/
addpath ~/new_matlab_functions/
initCobraToolbox;
%% calculate the original integration
[OFD_ori,~,~,~,~,~,~,~,~,~,sensitivity_ori,~,Nfit_latent_ori,minTotal_OFD_ori] = autoIntegration_latent(worm,true,storageCoef,sideCoef,epsilon_f,epsilon_r,capacity_f,capacity_r, ATPm, ExpCatag);
%% perturb the zero and 1 scores for latent sensitivity 
tol = 1e-7; % the default flux numerical tolerance is alway 1e-7
if clusterFlag
    % load cluster environment
    loadCluster(60,'1024','large','4:00')
else
    parpool(2)
    pctRunOnAll initCobraToolbox;
end
%make the output list
sensitivity_Name = sensitivity_ori.Name;
sensitivity_rxnID = sensitivity_ori.rxnID;
sensitivity_Score = sensitivity_ori.Score;
%make the rxn list 
ssInd = find(sensitivity_ori.Score <=1 & sensitivity_ori.Score >=-1);
sensitivity_Score_tmp = zeros(length(ssInd),1);

parfor idx = 1:numel(ssInd)
    i = ssInd(idx);
    model_ptb = worm;
    rxnInd = find(strcmp(model_ptb.rxns,sensitivity_rxnID{i}));
    if abs(OFD_ori(rxnInd)) < tol %zero -> non-zero; give a negative score
        if strcmp(sensitivity_Name{i}(end-1:end),'_f')
            % if this reaction is a forward reaction
            model_ptb.lb(rxnInd) = epsilon_f(rxnInd);
            try %in case infeasible
                [~,~,~,~,~,~,~,~,~,~,~,~,Nfit_latent_p,minTotal_OFD_p] = ...
                    autoIntegration_latent_manual(model_ptb,true,storageCoef,sideCoef,epsilon_f,epsilon_r,[],[], ATPm, ExpCatag,0,0,sensitivity_ori);%whether do minTotal depends on whther do full sensitivity analysis
                % scoring 
                if Nfit_latent_p < Nfit_latent_ori
                    sensitivity_Score_tmp(idx) = -1.5;
                elseif minTotal_OFD_p > minTotal_OFD_ori + 1e-7 %the numerical tolerance
                    sensitivity_Score_tmp(idx) = -1;
                else
                    sensitivity_Score_tmp(idx) = 0;
                end
            catch
                sensitivity_Score_tmp(idx) = -1.5; %assume infeasible as a 1.5 sensitivity
                InfLog = fopen('./InfCases_OFDss.txt','a+');%write the infeasible cases to a log file
                fprintf(InfLog,'The rxn %s is infeasible after perturbation; the original flux is %.2f\n',sensitivity_Name{i},OFD_ori(rxnInd));
                fclose(InfLog);
            end
        else %its a reverse reaction
            model_ptb.ub(rxnInd) = -epsilon_r(rxnInd);
            try
                [~,~,~,~,~,~,~,~,~,~,~,~,Nfit_latent_p,minTotal_OFD_p] = ...
                    autoIntegration_latent_manual(model_ptb,true,storageCoef,sideCoef,epsilon_f,epsilon_r,[],[], ATPm, ExpCatag,0,0,sensitivity_ori);%whether do minTotal depends on whther do full sensitivity analysis
                if Nfit_latent_p < Nfit_latent_ori
                    sensitivity_Score_tmp(idx) = -1.5;
                elseif minTotal_OFD_p > minTotal_OFD_ori + 1e-7 %the numerical tolerance
                    sensitivity_Score_tmp(idx) = -1;
                else
                    sensitivity_Score_tmp(idx) = 0;
                end        
            catch
                sensitivity_Score_tmp(idx) = -1.5;
                InfLog = fopen('./InfCases_OFDss.txt','a+');%write the infeasible cases to a log file
                fprintf(InfLog,'The rxn %s is infeasible after perturbation; the original flux is %.2f\n',sensitivity_Name{i},OFD_ori(rxnInd));
                fclose(InfLog);
            end
        end
    else %it is orginally non-zero flux
        if (strcmp(sensitivity_Name{i}(end-1:end),'_f') && OFD_ori(rxnInd) > 0) || (strcmp(sensitivity_Name{i}(end-1:end),'_r') && OFD_ori(rxnInd) < 0) %f reaction and ori flux is f || r reaction and ori flux is r
            model_ptb.lb(rxnInd) = 0;
            model_ptb.ub(rxnInd) = 0; %just put it zero
            try
                [~,~,~,~,~,~,~,~,~,~,~,~,Nfit_latent_p,minTotal_OFD_p] = ...
                    autoIntegration_latent_manual(model_ptb,true,storageCoef,sideCoef,epsilon_f,epsilon_r,[],[], ATPm, ExpCatag,0,0,sensitivity_ori);%whether do minTotal depends on whther do full sensitivity analysis
                if Nfit_latent_p < Nfit_latent_ori
                    sensitivity_Score_tmp(idx) = 1.5;
                elseif minTotal_OFD_p > minTotal_OFD_ori + 1e-7 %the numerical tolerance
                    sensitivity_Score_tmp(idx) = 1;
                else
                    sensitivity_Score_tmp(idx) = 0;
                end
            catch
                sensitivity_Score_tmp(idx) = 1.5; 
                InfLog = fopen('./InfCases_OFDss.txt','a+');%write the infeasible cases to a log file
                fprintf(InfLog,'The rxn %s is infeasible after perturbation; the original flux is %.2f\n',sensitivity_Name{i},OFD_ori(rxnInd));
                fclose(InfLog);
            end
        else %is the nonflux reverse reaction, put it to the corresponding epsilon flux
            if strcmp(sensitivity_Name{i}(end-1:end),'_f')
                model_ptb.lb(rxnInd) = epsilon_f(rxnInd);
            else
                model_ptb.ub(rxnInd) = -epsilon_r(rxnInd);
            end
            try %in case infeasible
                [~,~,~,~,~,~,~,~,~,~,~,~,Nfit_latent_p,minTotal_OFD_p] = ...
                    autoIntegration_latent_manual(model_ptb,true,storageCoef,sideCoef,epsilon_f,epsilon_r,[],[], ATPm, ExpCatag,0,0,sensitivity_ori);%whether do minTotal depends on whther do full sensitivity analysis
                % scoring 
                if Nfit_latent_p < Nfit_latent_ori
                    sensitivity_Score_tmp(idx) = -1.5;
                elseif minTotal_OFD_p > minTotal_OFD_ori + 1e-7 %the numerical tolerance
                    sensitivity_Score_tmp(idx) = -1;
                else
                    sensitivity_Score_tmp(idx) = 0;
                end
            catch
                sensitivity_Score_tmp(idx) = -1.5; %infeasible is also a strongest sensitivity
                InfLog = fopen('./InfCases_OFDss.txt','a+');%write the infeasible cases to a log file
                fprintf(InfLog,'The rxn %s is infeasible after perturbation; the original flux is %.2f\n',sensitivity_Name{i},OFD_ori(rxnInd));
                fclose(InfLog);
            end
        end
    end        
end
sensitivity_Score(ssInd) = sensitivity_Score_tmp;

sensitivity.Name = sensitivity_Name;
sensitivity.rxnID = sensitivity_rxnID;
sensitivity.Score = sensitivity_Score;
end
