% make the initial category of seq data
% the criterias are: 5% FDR, 2 fold change cutoff or no FC cutoff

%% part I: process the differential expression data
% because DE analysis doesnt apply to this dataset, so we change the
% strategy to picking the top 
% the original rowwise data is more normally distributed (more match CLT)
% (data not shown), while the rowwise log2(TPM+1) is more skewed. we use
% the z-transformation of the TPM matrix (rowwise Z) to determine if a gene
% is significantly differential expressed.
zCutoff = 2.33; % p < 0.01
load('TPM.mat')
expMat = cell2mat(TPM(2:end,2:end));
TPM = cell2table(TPM(2:end,:),'VariableNames',TPM(1,:));
zMat = normalize(expMat,2);
DEgenes = TPM.Gene(max(abs(zMat),[],2)>zCutoff);
% in the refinement, a gene with z score higher than cutoff will be upgrade
% to high if necessary;
save('DElist.mat','DEgenes','zMat');
%% part II: the classification based on absolute value
% only metabolic genes
fitData = TPM{:,2:end};
%fitData = [data.rep4_daf16;data.rep4_daf2;data.rep4_double;data.rep4_N2;data.rep5_daf16;data.rep5_daf2;data.rep5_double;data.rep5_N2;data.rep6_daf16;data.rep6_daf2;data.rep6_double;data.rep6_N2];
histogram(log2(fitData))
% fit bimodel guassian
rng(1126)
x = log2(fitData);
x(isinf(x)) = [];
fit = fitgmdist(x',2,'Options',statset('Display','final','MaxIter',3000,'TolFun',1e-9),'Replicates',10);
% note: the sigma in output is sigma^2
%% visualization
%options = statset('MaxIter',300, 'MaxFunEvals',600);
%paramEsts = mle(x, 'pdf',pdf_normmixture, 'start',start, ...
%                          'lower',lb, 'upper',ub, 'options',options)
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
%%
xline( fit.mu(2)-sqrt(fit.Sigma(2))*3,'--r');
xline(fit.mu(1)-sqrt(fit.Sigma(1))*3,'--r');
% xline(fit.mu(3) + sqrt(fit.Sigma(3)),'--r');
% xline(-1.71,'--r');
xline(min(x),'--r');

%% gene catagorization
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