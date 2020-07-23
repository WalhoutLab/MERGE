import pandas as pd
import json
from CatExp import *

##1-Loading the input.

#Files that have the following information are needed:
#   A table of gene expression data.
#   The set of genes in the metabolic network model.
#If the table and model use different gene annotations, then a map (dictionary) from table names to model names is also needed. 

#EXAMPLE - C. elegans
#Load expression table for seven tissues from Cao et al., 2017. Set Wormbase ID column (gene IDs) as the index.
Texp_cel=pd.read_csv('Input/CEL_Cao2017.tsv', sep='\t');
Texp_cel=Texp_cel.set_index('wbid');
#Load mapper that converts Wormbase IDs of iCEL1314 genes to gene names.
with open('Input/mapWBID2GeneName_iCEL1314.json') as json_file: 
     mapID2name_cel=json.load(json_file); 
#Get model gene list. In this case, it is just the values of the mapID2name_cel dictionary, since this dictionary covers only model genes and covers all of them.
icelGenes=mapID2name_cel.values();

#EXAMPLE - Human
#Load consensus expression table for 68 tissues from Human Protein Atlas (https://www.proteinatlas.org/about/download). Format the table such that rows are genes and columns are tissues.
Texp_hsa=pd.read_csv('Input/HSA_consensus.tsv', sep='\t');
Texp_hsa=Texp_hsa.pivot(index='Gene',columns='Tissue',values='NX');
#Load mapper that converts Ensemble IDs of human genes to HGNC IDs.
with open('Input/mapENS2HGNC.json') as json_file: 
     mapID2name_hsa=json.load(json_file); 
#Get gene list for Recon 2.2 model (provided).
with open('Input/ReconGenes.txt') as f:
    reconGenes = f.read().splitlines()



##2-Categorize model genes based on absolute cutoffs.

#In this step, each gene is going to be categorized as one of the following: 
#   High (highly expressed)
#   Moderate (modertaley expressed)
#   Low (lowly expressed)
#   Rare (rarely expressed).
#
#Expression level thresholds that define these categories will be determined based on the histogram of average gene expression. Two or three Gaussian curves will be fitted to this histogram and the mean and standard deviation of these curves will be used to define the thresholds.  

#EXAMPLE - C. elegans
#Get average expression as a column in the expression table.
Texp_cel['ave']=Texp_cel.mean(numeric_only=True, axis=1);
#Try fitting two Gaussian curves to the histogram of data (in log2 form, without zeros), using the bimodal() function. 
Stats_cel=bimodal(Texp_cel['ave'],expectedStats=[0,2.,5.,1.]); #Plots the histogram and fit on the screen.
#expectedStats is the expected mean and standard deviation of the first curve followed by those of the second. The algorithm uses these parameters as initial values during curve fitting. If the fit is not good seemingly because of the bad positioning of one of the curves, this step can be repeated by changing expectedStats with more educated guesses based on the plot.
#In this case, both visual inspection and Stats_cel['Rsquared'] statistic (>0.99) suggests that this is a good fit. Therefore, the means (mu) and standard deviations (sigma) of the two curves in Stats_cel can be used to determine thresholds in the next step.
#For example, the upper boundary for rarely expressed can be mu1, the upper boundary for lowly expressed mu1+sigma1, and the lower boundary for highly expressed mu2. 
#You can use the categorize_absCutoff() function with these three cutoffs in the respective order (as powers of 2, since the Stats_cel indicate the curve statistics in the logarithmic scale).
#But to avoid categorizing all genes in the expression data, first extract model genes from the expression table using the subTable() function.
Texp_cel_model=subTable(Texp_cel,icelGenes,mapID2name_cel);
#Then, use this table together with the abovementioned thresholds to obtain a table with genes categorized into the expression groups as defined above.
Tcat_icel=categorize_absCutoff(Texp_cel_model,2**np.array([Stats_cel['mu1'],Stats_cel['mu1']+Stats_cel['sigma1'],Stats_cel['mu2']]),excludeCols=['ave']);  #excludeCols are numerical columns to exclude from categorization
#As an example:
print Tcat_icel.loc[['algn-13']];
#gives the categories for the algn-13 gene:
#        Neurons     Gonad Hypodermis   Pharynx Body_wall_muscle      Glia Intestine
#algn-13     Low  Moderate        Low  Moderate         Moderate  Moderate      Rare

#EXAMPLE - Human
#See C. elegans example above for details.
#First, get average expression as a column in the expression table.
Texp_hsa['ave']=Texp_hsa.mean(numeric_only=True, axis=1);
#Try fitting two Gaussian curves to the histogram of data (in log2 form, without zeros), using the bimodal() function. 
Stats_hsa=bimodal(Texp_hsa['ave'],expectedStats=[0,2.,5.,1.]); #Plots the histogram and fit on the screen.
#In this case, both visual inspection and Stats_cel['Rsquared'] statistic (<0.99) suggests that this is not a good fit. It also seems that the large population has a shoulder. Thus, try trimodal fitting next.
Stats_hsa=trimodal(Texp_hsa['ave'],expectedStats=[-4,1,0,1,3.5,1]);
#where, expectedStats is the expected mean and standard deviation of the three fitted curves.
#The resulting fit is much better with R2~0.999 and visual inspection. The statistics of the fitted curves can be used for setting the thresholds (see below).
#Prepare the expression table of Recon 2.2 for categorization.
Texp_hsa_recon=subTable(Texp_hsa,reconGenes,mapID2name_hsa);
#Now, categorize the genes with the help of cutoffs from curve fitting stats (see C. elegans example above for the input to the categorize_absCutoff() function).
Tcat_recon=categorize_absCutoff(Texp_hsa_recon,2**np.array([Stats_hsa['mu1']-Stats_hsa['sigma1'],Stats_hsa['mu1'],Stats_hsa['mu3']]),excludeCols=['ave']);  #excludeCols are numerical columns to exclude from categorization



##3-Categorize moderately expressed genes based on relative expression levels.

#In this step, moderately expressed genes from the previous step are reevaluated based on relative expression levels and may be assigned to lowly or highly expressed categories. 
#For each gene, the relative expression algorithm first sorts tissues with respect to the expresison level of the gene (from low to high) and then evaluates the jumps, in fold change, from one tissue to the next.

#EXAMPLE - C. elegans
#Use relativeExp() function with expression value table and category table from the previous steps.
Tcat_icel_final=relativeExp(Texp_cel_model,Tcat_icel,2**(Stats_cel['mu2']-Stats_cel['sigma2']),fc_mid=1.5,fc_end=4.);
#where, the third term is yet another threshold that defines the lower boundary for assigning a gene to High category and also the upper boundary for assigning a gene to Low category. fc_mid is the minimum fold change value used to define a significant increase from one tissue to the next in the middle of the sorted expression profile (see above), and fc_end is the same except that it is applied to the increase from the first to the second and from the penultimate to the last tissue. 
#Tcat_icel_final describes the final categories of each gene in the model. As an example:
print Tcat_icel_final.loc[['algn-13']];
#gives the final categories for the algn-13 gene:
#        Neurons Gonad Hypodermis Pharynx Body_wall_muscle Glia Intestine
#algn-13     Low  High        Low    High              Low  Low      Rare
#
#To plot the ascending expresison profile of a with categorization, use plotCatExp function, which needs expression table without the average column:
Texp_cel_model_noave=Texp_cel_model.drop('ave',axis=1);
plotCatExp('E01A2.1',Texp_cel_model_noave,Tcat_icel,Tcat_icel_final);
#which plots the profile with color code: {'High':'green','Moderate':'gray','Low':'orange','Rare':'red'}, and hatch indicates categories that changed after relative expression analysis.
#Finally, to save the result:
Tcat_icel_final.to_csv('Output/Tcat_iCEL1314.tsv',sep='\t',index_label='Gene')

#EXAMPLE - Human
#Use relativeExp() function with expression value table and category table from the previous steps (see C. elegans example above for the description of input).
Tcat_recon_final=relativeExp(Texp_hsa_recon,Tcat_recon,2**(Stats_hsa['mu2']-Stats_hsa['sigma2']),fc_mid=2.,fc_end=4.);
#Tcat_recon_final describes the final categories of each gene in the model. 
#You can use the stackedCat() finction to see the categorical distribution of genes in each tissue.
stackedCat(Tcat_recon_final);
#To save the result
Tcat_recon_final.to_csv('Output/Tcat_recon.tsv',sep='\t',index_label='Gene');



