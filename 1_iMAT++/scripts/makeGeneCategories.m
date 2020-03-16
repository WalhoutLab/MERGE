%% Overview
% This is an interactive gene category builder that helps user to fit and
% build the category input for iMAT++. The expression quantifications (TPM
% or FPKM) could be provided as a csv table.
%% part I: load the expression data and adjust the format
TPM = readtable('./../inputs/TPM.csv');
% (specific to worm model)
% worm expression data are often labeled with WormBase ID (WBID). So, we
% provided a simple lookup table for iCEL1314.

% for other models or inputs, please make sure:
% 1. the columns are conditions to be analyzed
% 2. the rows are genes with gene ID used in the model
%% part II: the classification based on absolute value
% user can perform the fitting both by all genes or only metabolic genes
fitData = TPM{:,2:end};
histogram(log2(fitData))
% fit bimodel guassian
rng(1126)
x = log2(fitData);% this automatically ignored all the zeros
x(isinf(x)) = [];
fit = fitgmdist(x',2,'Options',statset('Display','final','MaxIter',3000,'TolFun',1e-9),'Replicates',10);
% note: the sigma in output is sigma^2
%% visualization of the guassian distribution fitted
bins = -15:.5:15;
h = bar(bins,histc(x,bins)/(length(x)*.5),'histc');
h.FaceColor = [.9 .9 .9];
xgrid = linspace(1.1*min(x),1.1*max(x),200);
pdfgrid = pdf(fit,xgrid');
%pdfgrid = pdf_normmixture(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
hold on
plot(xgrid,pdfgrid,'-')
hold off
xlabel('x')
ylabel('Probability Density')
%% label the desired cutoff
xline( fit.mu(2)-sqrt(fit.Sigma(2))*3,'--r');
xline(fit.mu(1)-sqrt(fit.Sigma(1))*3,'--r');
% xline(fit.mu(3) + sqrt(fit.Sigma(3)),'--r');
% xline(-1.71,'--r');
xline(min(x),'--r');

%% build the gene catagories
zero2low = min(x);%fit.mu(1)-sqrt(fit.Sigma(1))*3;%-sqrt(fit.Sigma(1))*2;%considering the absolute TPM is still high for mu1, we use u1-2sigma to define zero category
dynamic2high = fit.mu(2)-sqrt(fit.Sigma(2))*3;%very greedy for high, call all detectable met genes as high
names = TPM.Properties.VariableNames(2:end);
load('TPM.mat')
TPM = cell2table(TPM(2:end,:),'VariableNames',TPM(1,:));
metgenes = model.genes;
for myName = names
    myTPM = log2(TPM.(myName{:}));
    GeneID = TPM.Gene;
    ExpCatag.zero = GeneID(myTPM < zero2low);
    %ExpCatag.low = GeneID(TPM >= zero2low & TPM < low2dynamic);
    ExpCatag.dynamic = GeneID(myTPM >= zero2low & myTPM < dynamic2high);
    ExpCatag.high = GeneID(myTPM >= dynamic2high);
    % the uncalled genes are in dynamic
    ExpCatag.dynamic = [ExpCatag.dynamic; metgenes(~ismember(metgenes,GeneID))];
    save(['categ_',myName{:},'.mat'],'ExpCatag');
end
%% analyze if the thereshold makes sense
myTPM = log2(TPM.(names{1}));
GeneID = TPM.Gene;
lowinall = GeneID(myTPM < zero2low);
highinall = GeneID(myTPM >= dynamic2high);
highmerge = GeneID(myTPM >= dynamic2high);

for myName = names
    myTPM = log2(TPM.(myName{:}));
    GeneID = TPM.Gene;
    ExpCatag.zero = GeneID(myTPM < zero2low);
    %ExpCatag.low = GeneID(TPM >= zero2low & TPM < low2dynamic);
    ExpCatag.dynamic = GeneID(myTPM >= zero2low & myTPM < dynamic2high);
    ExpCatag.high = GeneID(myTPM >= dynamic2high);
    % the uncalled genes are in dynamic
    ExpCatag.dynamic = [ExpCatag.dynamic; metgenes(~ismember(metgenes,GeneID))];
    lowinall = intersect(lowinall,ExpCatag.zero);
    highinall = intersect(highinall,ExpCatag.high);
    highmerge = union(highmerge,ExpCatag.high);
end