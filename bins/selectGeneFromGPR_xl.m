function expressionCol = selectGeneFromGPR_xl(model, gene_names, gene_exp, parsedGPR)
% Map gene expression to reaction expression using the GPR rules. An AND
% will be replaced by max and an OR will be replaced by min.
%
% USAGE:
%   expressionCol = selectGeneFromGPR(model, gene_names, gene_exp, parsedGPR)
%
% INPUTS:
%   model:          COBRA model struct
%   gene_names:     gene identifiers corresponding to gene_exp. Names must
%                   be in the same format as model.genes (column vector)
%                   (as returned by "findUsedGeneLevels.m")
%   gene_exp:       quodra-level expression
%   parsedGPR:      GPR matrix as returned by "GPRparser.m"
%
% OUTPUTS:
%   expressionCol:  reaction expression, corresponding to model.rxns.
%                   No gene-expression data and orphan reactions will
%                   be given a value of -1.
%
% AUTHOR: Anne Richelle, May 2017
% Xuhang Li, 08192019 - optimized for speed 
expressionCol = -1*ones(length(model.rxns),1); %-1 means unknown/no data
    for i = 1:length(model.rxns)
        curExprArr=parsedGPR{i};
        curExpr= nan(length(curExprArr),1);
        for j=1:length(curExprArr)
            if length(curExprArr{j})>=1
                geneID = ismember(gene_names,curExprArr{j});
                if any(geneID) %if the gene is measured
                     curExpr(j)= min(gene_exp(geneID)); %If there is data for any gene in 'AND' rule, take the minimum value
                end
            end
        end
        if any(~isnan(curExpr))
            expressionCol(i)=max(curExpr);%If there is data for any gene in the 'OR' rule, take the maximum value
        end
    end
    %now, all the reactions are labeled. The Zero expression will overide
    %one in AND gate and be overided by one in OR and so on
    %next, we retrieve all the reactions associated with high expression
    %genes(level = 3), and only exluded low and zero reactions. Therefore, high
    %reaction will overide moderate level genes.
    HGenes = gene_names(gene_exp == 3);
    HrxnsInd = any(model.rxnGeneMat(:,ismember(model.genes,HGenes)),2);
    NonLowZero = expressionCol ~= 0 & expressionCol ~= 1; %all orphan, high labeled, moderate labeled
    HrxnsCandidate = HrxnsInd & NonLowZero;
    expressionCol(HrxnsCandidate) = 3; %all HGene associated reactions (no conflict with low and zero) are re-labeled as level = 3
end



