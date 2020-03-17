This guidance shows how to reproduce the flux prediction in Fig 2D by the original iMAT algorithm from [Shlomi et al., 2008](https://pubmed.ncbi.nlm.nih.gov/18711341/) and as modified in [Yilmaz and Walhout, 2016](https://pubmed.ncbi.nlm.nih.gov/27211857/) with a flux minimization in the end. 

### Running the integration for the dual-tissue model

This script helps run iMAT for the seven <i>C. elegans</i> tissues for which aggregated gene expression data is available from [Cao <i>et al</i>., 2017](https://pubmed.ncbi.nlm.nih.gov/28818938/).Input files used are indicated in the walkthrough script [here](myFlux.m). In short, these include [the dual version of iCEL1314](./../input/Tissue.mat), [premade gene categories](./../input/geneCategories.json) and [premade epsilon values](./../input/epsilon.json). 

For converting the excel-formatted model (<i>i.e.</i>, the table downloaded from supplemental material) to CORRA-version model, see [here](makeWormModel.m). For categorizing gene expression data as highly, lowly, rarely, and moderately expressed genes, please see [our paper](paper link) that describes the categorization algorithm in detail. The categories used in this study is provided in json format, and can be loaded into MATLAB following instructions in the walkthrough script. Users can view the categories after loading the file. To obtain a rough gene category from any expression profile, please follow [the gene category generator](./scripts/makeGeneCategories.m). Note that the output from this generator is different from the fine-tuned category in this study, because of lack of the heuristic refinement. Gene expression can be categorized using other methods as well (<i>e.g.</i>, with thresholds on percentiles, such as selecting top 10% percentile as highly expressed genes), depending on the dataset and user preferences. 

Finally, for defining the epsilon values (flux thresholds for each reaction that discriminates a real flux from zero flux during the integration), please refer to [epsilon generator](./../bins/makeEpsilonSeq.m). The method for defining epsilons can also vary from user to user. For example, a standard epsilon value of 0.001 may be used for all reactions. But we do not recommend this approach as some reactions tend to carry low flux and forcing values of 0.001 on these reactions cause excessive flux values (>1000) in other parts of the network due to mass balance rules.

To reproduce the flux distribution prediction, simply run:
```
matlab < myFlux.m
```
Output (Flux distribution and other metrics) will be written in output/TissueName.mat

Please refer to the [integration function](iMAT_xl.m) for the output format.

For more instructions of running the pipeline in a custom dataset, please refer to the header of each function included in the walkthrough script.

### Important Note on Computational Performance

If you are intended to run iMAT algorithm to reproduce the results from this study, please make sure the programs are excuted in a high-performance workstation (we suggest >= 20 CPU cores and > 20 Gb memory). In our experience, though iMAT++ is fast enough to be run in a laptop even with a large network model such as the dual-tissue model, the original iMAT is computationally intensive. It is typical to take 30mins to hours to complete the integration for all tissues in a 20-core workstation. 

