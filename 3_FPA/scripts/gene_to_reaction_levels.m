function [reaction_levels, reaction_status] = gene_to_reaction_levels( model, genes, levels, f_and, f_or )
% Convert gene expression levels to reaction levels using GPR associations.
% Level is NaN if there is no GPR for the reaction or no measured genes.
%
% INPUTS
%       model - cobra model
%       genes - gene names
%       levels - gene expression levels
%       f_and - function to replace AND
%       f_or - function to replace OR
%
% OUTPUTS
%       reaction_levels - reaction expression levels
%
% Author: Daniel Machado, 2013
% adapted. the output is a string of expression levels for logic AND blocks

    reaction_levels = cell(length(model.rxns), 1);
    reaction_status = 1000*ones(length(model.rxns), 1);
    preParsedGrRules = preparseGPR_xl(model.grRules);
    for i = 1:length(model.rxns)
        [level,status] = eval_gpr(preParsedGrRules{i}, genes, levels, f_and, f_or);
        reaction_levels{i} = level;
        reaction_status(i) = status;
    end
end

