%% Overview
% This is an interactive gene category builder that helps user to fit and
% build the categories for iMAT++. The expression quantifications (TPM
% or FPKM) could be provided as a csv table.
% ..Author: Xuhang Li, Mar 2020
%% part I: load the expression data and adjust the format
% we use the RNA-seq data from Bulcha et al, Cell Rep (2019) as an example
TPM = readtable('./input/exampleExpression.csv');
% (specific to worm model)
% worm expression data are often labeled with WormBase ID (WBID). So, we
% provided a simple lookup table for iCEL1314.
lookupTbl = readtable('./input/IDtbl.csv');
metGenesInd = ismember(TPM.GeneID,lookupTbl.WormBase_Gene_ID);
wormbaseTbl = readtable('./input/wormbasetable.tsv','FileType','text');
proteinCodingGenes = wormbaseTbl.WormBaseID(strcmp(wormbaseTbl.Type,'protein coding'));
proteinCodingInd = ismember(TPM.GeneID,proteinCodingGenes);
% for other models or inputs, please make sure:
% 1. the columns are conditions to be analyzed
% 2. the rows are genes with gene ID used in the model
%% part II: the classification based on absolute value
% Rationale:
% user can perform the fitting both by all protein coding genes or only
% metabolic genes.

% Please note that the distribution of TPM (FPKM) varies from dataset to
% dataset, so it is not necessary that the categorization thresholds we used in
% the tissue modeling is uniformly best choices for all dataset. In fact, our
% experience suggests that it is best to use ALL PROTEIN CODING GENES for
% SINGLE CELL RNA-SEQ DATA, but just METABOLIC GENES for BULK RNA-SEQ data,
% in C. elegans. The possible reason is that the bulk RNA-seq of the whole
% animal lacks the resolution to seperate the lowly expressed, tissue-specific genes
% from the expression noise (the actual "lowly expressed subpop") into distinct subpopulations.
% However, because metabolic genes are generally highly expressed, so the actual
% "lowly expressed subpop" still gets seperated if only uses metabolic genes. 

% fitData = TPM{proteinCodingInd,5:end};
fitData = TPM{metGenesInd,5:end};% we perform guassian fitting only on metabolic genes
histogram(log2(fitData))
% fit bimodel guassian
rng(1126)
x = log2(fitData);% this automatically ignored all the zeros
x(isinf(x)) = [];
fit = fitgmdist(x',2,'Options',statset('Display','final','MaxIter',3000,'TolFun',1e-9),'Replicates',10);
% note: the sigma in output is sigma^2
%% visualization of the guassian distribution fitted
% user can inspect how well the fitting is and whether it seperates into
% two subpopulation as expected.
bins = -15:.5:15;
h = bar(bins,histc(x,bins)/(length(x)*.5),'histc');
h.FaceColor = [.9 .9 .9];
xgrid = linspace(1.1*min(x),1.1*max(x),200);
pdfgrid = pdf(fit,xgrid');
hold on
plot(xgrid,pdfgrid,'-')
hold off
xlabel('x')
ylabel('Probability Density')
%% label the desired cutoff
xline( fit.mu(2),'--r');
xline(fit.mu(1),'--r');
xline(fit.mu(2) + sqrt(fit.Sigma(2)),'--r');
% though we find that the mean of the two subpopulation gives best
% categories in the datasets we tested, we recommend users to interactively
% evaluate their dataset and find the best thresholds
%% build the gene catagories
TPM = readtable('./input/exampleExpression.csv');
zero2low = fit.mu(2);%set thresholds
low2dynamic = fit.mu(2) + sqrt(fit.Sigma(2));%set thresholds
dynamic2high = fit.mu(1);%set thresholds
names = TPM.Properties.VariableNames(5:end);%choose all the sample names
metgenes = model.genes;
TPM = TPM(ismember(TPM.GeneID,lookupTbl.WormBase_Gene_ID),:);%we only keep the genes in the model
 % specifically, we need to convert the gene ID for C. elegans model
[A,B] = ismember(TPM.GeneID,lookupTbl.WormBase_Gene_ID);
TPM.GeneID(A) = lookupTbl.ICELgene(B(A));
for myName = names
    myTPM = log2(TPM.(myName{:}));
    GeneID = TPM.GeneID;
    ExpCateg.zero = GeneID(myTPM < zero2low);
    ExpCateg.low = GeneID(myTPM >= zero2low & myTPM < low2dynamic);
    ExpCateg.dynamic = GeneID(myTPM >= low2dynamic & myTPM < dynamic2high);
    ExpCateg.high = GeneID(myTPM >= dynamic2high);
    % the uncalled genes (i.e., NA and ND)are in dynamic (moderately expressed)
    ExpCateg.dynamic = [ExpCateg.dynamic; metgenes(~ismember(metgenes,GeneID))];
    save(['input/exampleGeneCategories/categ_',myName{:},'.mat'],'ExpCateg');
end
%% analyze if the thereshold makes sense
% A critical signature of good gene categories is that most metabolic genes are
% called as high when merging the high categories for all conditions
% (assuming the metabolism across conditions varies), while only a core set
% of genes are high in all of the conditions. This indicates that the
% categorization captures the expression regulation of metabolic genes
% across conditions. Similarly, we expect the zero category to have the
% similar feature. 
myTPM = log2(TPM.(names{1}));
GeneID = TPM.GeneID;
ZeroInAll = GeneID(myTPM < zero2low);
ZeroMerge = GeneID(myTPM < zero2low);
HighInAll = GeneID(myTPM >= dynamic2high);
HighMerge = GeneID(myTPM >= dynamic2high);
N_zero = [];
N_low = [];
N_dynamic = [];
N_high = [];
for myName = names
    myTPM = log2(TPM.(myName{:}));
    GeneID = TPM.GeneID;
    ExpCateg.zero = GeneID(myTPM < zero2low);
    ExpCateg.low = GeneID(myTPM >= zero2low & myTPM < low2dynamic);
    ExpCateg.dynamic = GeneID(myTPM >= low2dynamic & myTPM < dynamic2high);
    ExpCateg.high = GeneID(myTPM >= dynamic2high);
    % the uncalled genes (i.e., NA and ND)are in dynamic (moderately expressed)
    ExpCateg.dynamic = [ExpCateg.dynamic; metgenes(~ismember(metgenes,GeneID))];
    ZeroInAll = intersect(ZeroInAll,ExpCateg.zero);
    ZeroMerge = union(ZeroMerge,ExpCateg.zero);
    HighInAll = intersect(HighInAll,ExpCateg.high);
    HighMerge = union(HighMerge,ExpCateg.high);
    N_zero = [N_zero;length(ExpCateg.zero)];
    N_low = [N_low;length(ExpCateg.low)];
    N_dynamic = [N_dynamic;length(ExpCateg.dynamic)];
    N_high = [N_high;length(ExpCateg.high)];
end
fprintf('%d/%d are highly expressed genes in all conditions\n',length(HighInAll),length(model.genes));
fprintf('%d/%d are highly expressed genes in at least one condition\n',length(HighMerge),length(model.genes));
fprintf('%d/%d are rarely expressed genes in all conditions\n',length(ZeroInAll),length(model.genes));
fprintf('%d/%d are rarely expressed genes in at least one conditions\n',length(ZeroMerge),length(model.genes));
% we offer an additional QC figure for category making (similar to fig. 1E
figure(3)
stackN = [N_zero,N_low,N_dynamic,N_high];
bar(1:60,stackN,'stacked')
xlabel('cell line No.')
legend({'zero','low','dynamic','high'});
%% further considerations
% Now you already have a rough gene category for each conditions. However,
% as mentioned in the paper, we further used a heuristic algorithm to
% refine the moderately expressed genes (aka, dynamic category).
% However, the heuristic algoritm works best with our tissue modeling task 
% as we use single cell sequencing data. For general application on bulk RNA-seq 
% dataset, we recommand to start with the rough category, or refine the category based on
% pair-wise differential expression (user need to develop their own code for this).