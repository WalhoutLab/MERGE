This guidance shows how to reproduce the flux potential analysis (FPA) result in Fig X by FPA. Additionally, we provided a walkthrough guidence for running FPA on the generic model of iCEL1314 as well as any other metabolic network model.

### Run the integration for dual model

This script helps run FPA for all seven C. elegans tissues. Used input files are indicated in the walkthrough script [here](myFlux.m). In short, [the dual version of iCEL1314](./../input/Tissue.mat), [premade gene categories](./../input/geneCategories.json) and [premade epsilon values](epsilon.json) are used. For converting the excel-format model to CORRA-version model, see [here](makeWormModel.m); For making gene categories, please see [our paper](paper link) for detailed method to obtain a fine-tuned gene category. The categories used in this study is provided in json format, and can be loaded into MATLAB following instructions in the walkthrough script. Users can view the categories after loading the file. To obtain a rough gene category from any expression profile, please follow [the gene category generator](To be uploaded). Note that the output from this generator is different from the fine-tuned category in this study, because of lack of the heuristic refinement. Finally, for making the epsilon sequencing, please refer to [epsilon generator](./../bins/makeEpsilonSeq.m).

To reproduce the flux distribution prediction, simply run:
```
matlab < writeFluxPotential_dual.m
```
Output(Flux distribution and other metrics) will be written in output/TissueName.mat

Please refer to the [integration function](scripts/autoIntegration_latent.m) for output format.

For more instructions of running the pipeline in a custom dataset, please refer to the header of each functions included in the walkthrough script.
