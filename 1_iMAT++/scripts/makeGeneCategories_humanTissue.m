%% Overview
% This is an interactive gene category builder that helps user to fit and
% build the category input for iMAT++. The expression quantifications (TPM
% or FPKM) could be provided as a csv table.
% ..Author: Xuhang Li, Mar 2020
%% part I: load the expression data and adjust the format
% we use the "RNA HPA tissue gene data" for human tissues downloaded from https://www.proteinatlas.org/about/download
TPM = readtable('./input/humanModel/rna_tissue_hpa.xlsx');% TPMs for all the protein coding genes provided in a processed table
% the format adjustment procedure varies from data to data, generally we need to make sure:
% 1. the columns are conditions to be analyzed
% 2. the rows are genes with gene ID used in the model (gene label should be unique)

% the following is the codes we used to adapt human tissue RNA-seq data
% generate the expression matrix for all genes in the model 
load('./input/humanModel/Recon2_2.mat');% load model
GOI = intersect(model.genes,TPM.HGNC);
Tissues = unique(TPM.Tissue);
Gene = GOI;
output=table(Gene);
tmp = array2table(nan(length(Gene),length(Tissues)));
varNames = regexprep(Tissues,' |,','_');
varNames = regexprep(varNames,'-','_');
tmp.Properties.VariableNames = varNames';
output = [output,tmp];
for j = 2:size(output,2)
    tissueInd = strcmp(TPM.Tissue,Tissues{j-1});
    geneLabel = TPM.HGNC(tissueInd);
    TPMval = TPM.TPM(tissueInd);
    [A B] = ismember(output.Gene,geneLabel);
    output{A,j} = TPMval(B(A));
end
expMat = output{:,2:end};

% make the same matrix for all the genes
Gene_all = unique(TPM.HGNC);
output_all=table(Gene_all);
tmp = array2table(nan(length(Gene_all),length(Tissues)));
tmp.Properties.VariableNames = varNames';
output_all = [output_all,tmp];
for j = 2:size(output_all,2)
    tissueInd = strcmp(TPM.Tissue,Tissues{j-1});
    geneLabel = TPM.HGNC(tissueInd);
    TPMval = TPM.TPM(tissueInd);
    [A B] = ismember(output_all.Gene_all,geneLabel);
    output_all{A,j} = TPMval(B(A));
end
expMat_all = output_all{:,2:end};
%% part II: the classification based on absolute value
% it is very important to notice that the distribution of TPM (FPKM) varies from dataset to
% dataset, so it is not necessary the categorization thresholds we used in
% the tissue modeling is uniformly best choices for all dataset. User needs
% to inspect their own data for making the optimal decision. Indeed, in the
% human tissue demo, change the way of setting cutoffs (see below).

figure(1)
hold on
histogram(log2(mean(expMat,2)))% only metabolic genes
histogram(log2(mean(expMat_all,2)))% all protein coding genes

% fit bimodel guassian
fitData = mean(expMat_all,2);
rng(1126)
x = log2(fitData);% this automatically ignored all the zeros
x(isinf(x)) = [];
fit = fitgmdist(x,2,'Options',statset('Display','final','MaxIter',3000,'TolFun',1e-9),'Replicates',10); % we found three subpopulation gives best fit
% note: the sigma in output is sigma^2
%% visualization of the guassian distribution fitted
% user can inspect how well the fitting is
figure(2)
bins = -5:.5:12;
h = bar(bins,histc(x,bins)/(length(x)*.5),'histc');
h.FaceColor = [.9 .9 .9];
xgrid = linspace(1.1*min(x),1.1*max(x),200);
pdfgrid = pdf(fit,xgrid');
hold on
plot(xgrid,pdfgrid,'-')
hold off
xlabel('x')
ylabel('Probability Density')
%% label the desired cutoff and visualize
% this set of cutoffs give most meaningful categories in the following QC
xline(fit.mu(1),'--r');% dynamic/high
xline(fit.mu(1) - 2*sqrt(fit.Sigma(1)),'--r');% low/dynamic
xline(fit.mu(2) - 2*sqrt(fit.Sigma(2)),'--r');% zero/low

figure(1)
hold on
histogram(log2(mean(expMat,2)))
xline(fit.mu(2) - 2*sqrt(fit.Sigma(2)),'--r');
xline(fit.mu(1) - 2*sqrt(fit.Sigma(1)),'--r');
xline(fit.mu(1),'--r');
%% build the gene catagories
zero2low = fit.mu(2) - 2*sqrt(fit.Sigma(2));% set thresholds
low2dynamic =  fit.mu(1) - 2*sqrt(fit.Sigma(1));% set thresholds
dynamic2high = fit.mu(1);% set thresholds
names = output.Properties.VariableNames(2:end);% choose all the sample names
metgenes = model.genes;
for myName = names
    myTPM = log2(output.(myName{:}));% we only categorize the metabolic genes
    GeneID = output.Gene;
    ExpCateg.zero = GeneID(myTPM < zero2low);
    ExpCateg.low = GeneID(myTPM >= zero2low & myTPM < low2dynamic);
    ExpCateg.dynamic = GeneID(myTPM >= low2dynamic & myTPM < dynamic2high);
    ExpCateg.high = GeneID(myTPM >= dynamic2high);
    % the uncalled genes (i.e., NA and ND)are in dynamic (moderately expressed)
    ExpCateg.dynamic = [ExpCateg.dynamic; metgenes(~ismember(metgenes,GeneID))];
    save(['input/humanModel/categ_',myName{:},'.mat'],'ExpCateg');
end
%% analyze if the thereshold makes sense
% A critical signature of good gene categories is that most metabolic genes are
% called as high when merging the high categories for all conditions
% (assuming the metabolism across conditions varies), while only a core set
% of genes are high in all of the conditions. This indicates that the
% categorization captures the expression regulation of metabolic genes
% across conditions. Similarly, we expect the zero category to have the
% same feature. 
myTPM = log2(output.(names{1}));
GeneID = output.Gene;
ZeroInAll = GeneID(myTPM < zero2low);
ZeroMerge = GeneID(myTPM < zero2low);
HighInAll = GeneID(myTPM >= dynamic2high);
HighMerge = GeneID(myTPM >= dynamic2high);
N_zero = [];
N_low = [];
N_dynamic = [];
N_high = [];
for myName = names
    myTPM = log2(output.(myName{:}));
    GeneID = output.Gene;
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
% we offer an additional QC figure for category making
figure(3)
stackN = [N_zero,N_low,N_dynamic,N_high];
bar(1:length(Tissues),stackN,'stacked')
xlabel('cell line No.')
legend({'zero','low','dynamic','high'});
%% further considerations
% Now you already have a rough gene category for each conditions. However,
% as mentioned in the paper, we further used a heuristic algorithm to
% refine the moderately expressed genes (aka, dynamic category).
% The limitation for it is that the heuristic algoritm works best with our
% tissue modeling task as we use single cell sequencing data. For general
% application on bulk RNA-seq dataset, we recommend to start with the rough
% category, or refine the category based on pair-wise differential expression.