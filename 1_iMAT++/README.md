This guidance shows how to reproduce the flux prediction in Fig 3A by iMAT++. Additionally, we provided user-friendly walkthrough guidance to run iMAT++ with any given dataset or metabolic model.

### Reproducing the integration for C. elegans' dual tissue model

This script helps run iMAT++ for all seven C. elegans tissues. Used input files are indicated in the walkthrough script [here](myFlux.m). In short, [the dual version of iCEL1314](./../input/Tissue.mat), [premade gene categories](./../input/geneCategories.json) and [premade epsilon values](epsilon.json) are used. 
For converting the excel-format model (i.e., the table downloaded from supplemental material) to CORRA-version model, see [here](makeWormModel.m); For making gene categories, please see [our paper](paper link) for detailed method to obtain a fine-tuned gene category. The categories used in this study is provided in json format, and can be loaded into MATLAB following instructions in the walkthrough script. Users can view the categories after loading the file. To obtain a rough gene category from any expression profile, please follow [the gene category generator](./scripts/makeGeneCategories.m). Note that the output from this generator is different from the fine-tuned category in this study, because of lack of the heuristic refinement. Finally, for making the epsilon sequencing, please refer to [epsilon generator](./../bins/makeEpsilonSeq.m) (guidance provided in the walkthrough script for general application of iMAT++).

To reproduce the flux distribution prediction, simply run:
```
matlab < myFlux.m
```
Output(Flux distribution and other metrics) will be written in output/TissueName.mat

Please refer to the [integration function](scripts/autoIntegration_latent.m) for output format.

### Run iMAT++ analysis on generic C. elegans model with a given expression dataset

For more instructions of running the pipeline in a custom dataset, please refer to the header of each functions included in the walkthrough script.

### Run iMAT++ analysis on other metabolic model with a given expression dataset

We provided a walkthrough guidance of how to run the iMAT++ analysis on human metabolic model RECON 2.2. The expression profile of NCI-60 Cancer Cell lines [Reinhold WC et all, 2019](https://cancerres.aacrjournals.org/content/79/13/3514.long) was used as an example.
/*To be completed.*/
