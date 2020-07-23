
This guidance shows how to use CatExp.py tool to categorize genes in a metabolic network model. 

Starting with a table of all detected genes and their expression elevels in different tissues (or coditions), the final outcome is a table of metabolic genes categorized into four levels of expression: Rare, Low, Moderate, and High.

Step-by-step instructions are provided in walkthrough.py using example categorizations with <i>C. elegans</i> and human datasets and models. Inputs are provided in [Input](./Input/) folder:
  - CEL_Cao2017.tsv: A table of gene expression values in seven tissues of <i>C. elegans</i> from Cao <i>et al.</i> (2017).
  - mapWBID2GeneName_iCEL1314.json: A dictionary that maps Wormbase IDs of all genes in iCEL1314 metabolic network model of <i>C. elegans</i> to gene names used in this model.
  - HSA_consensus.tsv: A table of consensus gene expression values in 62 human tissues downloaded from https://www.proteinatlas.org/humanproteome/tissue. 
  - mapENS2HGNC.json: A dictionary that maps Ensemble IDs of human genes to HGNC IDs.
  - ReconGenes.txt: A list of genes used in human model Recon 2.2.
  
Output folder includes the results (<i>i.e.</i>,tables of categorized genes) from <i>C. elegans</i> and human examples.  


### Categorization based on absolute values of expression

The function categorize_absCutoff() is used to categorize genes based on three expression level cutoffs: rare cutoff, low cutoff, and high cutoff.     

Genes with expression levels less than rare cutoff are categorized as Rare. 
Other genes with expression levels less than low cutoff are categorized as Low.
Genes with expression levels greater than high cutoff are categorized as High.
All other genes are categorized as Moderate.

These cutoffs can be determined using gene expression histograms. Arbitrary cutoffs (<i>e.g.</i>, 90th percentile for high cutoff) can be used or the histogram can be statistically analyzed to determine the thresholds according to the peaks, shoulders, and tails observed. 

For the statistical analysis of histograms, CatExp provides two curve fitting tools that deconvolute subpopulations of genes with different log-normal distributions. The first is bimodal() function that fits two Gaussian curves, as superimposed, to the histogram of logarithm (base 2) of average gene expression levels across all tissues (conditions). Zero values are excluded from the fitting. The mean and standard deviations of the two subpopulations can then be used to determine cutoffs (<i>e.g.</i>, the mean of the high expression subpopulation can be used as the high cutoff). The second tool is trimodel() function, which does the same job as bimodel(), but with three curves representing three log-normally distributed subpopulations, instead of two curves.


### Categorization based on relative expression

The function relativeExp() is used to categorize genes based on their relative expression levels across the tissues (conditions). The absolute categorizations from the first step (see above) is an input for this function. In the end, some of the Moderate genes are recategorized by a heuristic algorithm as High or Low, depending on relative expression levels. See the help text of relativeExp() for details.    


### Plotting

CatExp also provides three plotting functions for the visualization of data.
  - Histogram plotters embedded in bimodel() and trimodel() functions.
  - stackedCat() function, which can take a table of categories (genes as rows and tissues or conditions as columns) and plots a stacked bar graph of categories (Rare, Low, Moderate, High) for each column.
  - plotCatExp() function, which receives a particular gene name, the expression level table, the categorization table from absolute values analysis, and the categorization table from relative values analysis to plot a bar chart that shows ascending tissue (or condition) profile of the gene based on expression level table. Categories are color coded and categories that changed during relative expresison analysis are highlighted.


