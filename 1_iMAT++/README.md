This guidance shows how to reproduce the flux prediction in Figures 2C, 2D, and 3A by iMAT++. Additionally, we provided user-friendly walkthrough guidance to run iMAT++ with any given dataset or metabolic model.

### Reproducing the integration for C. elegans' dual tissue model

This script helps run iMAT++ for seven <i>C. elegans</i> tissues for which aggregated gene expression data is available from [Cao <i>et al</i>., 2017](https://pubmed.ncbi.nlm.nih.gov/28818938/). Input files used are indicated in the walkthrough script [here](myFlux.m). In short, these include [the dual version of iCEL1314](./../input/Tissue.mat), [premade gene categories](./../input/geneCategories.json) and [premade epsilon values](./../input/epsilon.json). 

For converting the excel-formatted model (<i>i.e.</i>, the table downloaded from supplemental material) to CORRA-version model, see [here](makeWormModel.m). For categorizing gene expression data as highly, lowly, rarely, and moderately expressed genes, please see [our paper](paper link) that describes the categorization algorithm in detail. The categories used in this study is provided in json format, and can be loaded into MATLAB following instructions in the walkthrough script. Users can view the categories after loading the file. To obtain a rough gene category from any expression profile, please follow [the gene category generator](./scripts/makeGeneCategories.m). Note that the output from this generator is different from the fine-tuned category in this study, because of lack of the heuristic refinement. Gene expression can be categorized using other methods as well (<i>e.g.</i>, with thresholds on percentiles, such as selecting top 10% percentile as highly expressed genes), depending on the dataset and user preferences. Finally, for defining the epsilon values (flux thresholds for each reaction that discriminates a real flux from zero flux during the integration), please refer to [epsilon generator](./../bins/makeEpsilonSeq.m). 

To reproduce the flux distribution prediction, simply run:
```
matlab < myFlux.m
```
Output(Flux distribution and other metrics) will be written in output/TissueName.mat

Please refer to the [integration function](scripts/IMATplusplus.m) for output format.

### Run iMAT++ analysis on generic C. elegans model with a given expression dataset

We provided a step-by-step walkthrough script for running iMAT++ on the generic iCEL1314 with any given expression dataset. Users can start with this script to perform iMAT++ analysis on their own expression data. 

We recommand to go over this script interactively in MATLAB graphic interface. Please launch a MATLAB instance and set your current directory in this folder (1_iMAT++). Then, open and go through [the walkthrough script](walkthrough_generic.m) according to the instructions in it.

### Run iMAT++ analysis on other metabolic model with a given expression dataset

We provided a walkthrough guidance of how to run the iMAT++ analysis on human metabolic model RECON 2.2. The expression profile of NCI-60 Cancer Cell lines [Reinhold WC et all, 2019](https://cancerres.aacrjournals.org/content/79/13/3514.long) was used as an example.
/*To be completed.*/
