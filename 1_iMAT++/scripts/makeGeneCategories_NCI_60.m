%% Overview
% This is an interactive gene category builder that helps user to fit and
% build the category input for iMAT++. The expression quantifications (TPM
% or FPKM) could be provided as a csv table.
% ..Author: Xuhang Li, Mar 2020
%% part I: load the expression data and adjust the format
% we use the RNA-seq data from Reinhold WC et al., 2019 for NCI-60 RNA-seq as an example
FPKM = readtable('./input/humanModel/log2FPKM.csv');% this is the log2(FPKM+1)
%though it may not be necessary, we convert FPKM to TPM for cross-sample comparasion
TPM = FPKM;
TPM{:,3:end} = 2.^TPM{:,3:end}-1;
TPM{:,3:end} = FPKM{:,3:end} ./ sum(FPKM{:,3:end},1) * 1e6;
% the format adjustment procedure varies from data to data, generally we need to make sure:
% 1. the columns are conditions to be analyzed
% 2. the rows are genes with gene ID used in the model (gene label should be
% unique)
% the following is the codes we used to adapt raw RNA-seq data input
%generate the mat file 
load('./input/humanModel/Recon2_2.mat');%load model
GOI = intersect(model.genes,TPM.GeneID);
output(1,:) = [{'Gene'},TPM.Properties.VariableNames(3:end)];
for i = 1:length(GOI)
    output(i+1,1) = GOI(i);
    for j = 2:length(output(1,:))
        ind = ismember(TPM.GeneID,GOI(i));
        output(i+1,j) = {TPM.(output{1,j})(ind)};
        if sum(ind) ~= 1 %skip these non-uniquely mapped genes
            output(i+1,1) = {'Delete'};
        end
    end
end
output(strcmp(output(:,1),'Delete'),:) = [];
expMat = cell2mat(output(2:end,2:end));
output = cell2table(output(2:end,:),'VariableNames',output(1,:));
%% part II: the classification based on absolute value
% it is very important to notice that the distribution of TPM (FPKM) varies from dataset to
% dataset, so it is not necessary the categorization thresholds we used in
% the tissue modeling is uniformly best choices for all dataset. User needs
% to inspect their own data for making the optimal decision. Indeed, we
% change the thresholding criteria for avoiding falsely putting expressed
% genes into zero category

% visualize the distribution of all genes and metabolic genes
histogram(log2(mean(TPM{:,3:end},2)))% average for all genes
hold on
histogram(log2(mean(expMat,2)))% only metabolic genes
% fit bimodel guassian
fitData = mean(TPM{:,3:end},2);% fit by all genes
rng(1126)
x = log2(fitData);% this automatically ignored all the zeros
x(isinf(x)) = [];
fit = fitgmdist(x,2,'Options',statset('Display','final','MaxIter',3000,'TolFun',1e-9),'Replicates',10);
% note: the sigma in output is sigma^2
%% visualization of the guassian distribution fitted
% user can inspect how well the fitting did and whether it seperates into
% two subpopulation as expected.
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
%% label the desired cutoff
xline( fit.mu(2) -  2*sqrt(fit.Sigma(2)),'--r');%low/zero
xline(fit.mu(1),'--r');%dynamic/high
xline(fit.mu(2),'--r');%low/dynamic

figure(1)
histogram(log2(mean(TPM{:,3:end},2)))% all genes
hold on
histogram(log2(mean(expMat,2)))% only metabolic genes
xline(fit.mu(2) -  2*sqrt(fit.Sigma(2)),'--r');
xline(fit.mu(1),'--r');
xline(fit.mu(2),'--r');
% notice that the mu(2) is still a large number (~2^2.5 = 6) which should
% not be representing unexpressed genes. So we altered the thresholding
% criteria for zero category
%% build the gene catagories
zero2low = fit.mu(2) - 2*sqrt(fit.Sigma(2));%set thresholds
low2dynamic = fit.mu(2);%set thresholds
dynamic2high = fit.mu(1);%set thresholds
names = output.Properties.VariableNames(2:end);%choose all the sample names
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
bar(1:60,stackN,'stacked')
xlabel('cell line No.')
legend({'zero','low','dynamic','high'});
%% further considerations
% Now you already have a rough gene category for each conditions. However,
% as mentioned in the paper, we further used a heuristic algorithm to
% further refine the moderately expressed genes (aka, dynamic category).
% The limitation for it is that the heuristic algoritm works best with our
% tissue modeling task as we use single cell sequencing data. For general
% application on bulk RNA-seq dataset, we recommand to start with the rough
% category, or use the "refinementTool.m" to refine the category based on
% pair-wise differential expression.