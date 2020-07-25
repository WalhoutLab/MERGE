# Modeling tissue-relevant <i>C. elegans</i> metabolism at network, pathway, reaction, and metabolite levels
We developped the package, MERGE (MEtabolic models Reconciled with Gene Expression), to systematically integrate tissue-level gene expression data from single-cell RNAseq provided by [Cao <i>et al</i>., 2017](https://pubmed.ncbi.nlm.nih.gov/28818938/) with the metabolic network model of <i>C. elegans</i> named [iCEL1314](http://wormflux.umassmed.edu/index.html). MERGE is applicable to other expression datasets and other metabolic network models. 

The package includes three main programs: (1) matlab implementation of [iMAT++ algorithm](1_iMAT++), which semi-quantitatively integrates categorized gene expression data with the model at system-level; (2) matlab implementation of [original iMAT algorithm](2_iMAT) (for comparison and reproduction of results in the paper), and (3) matlab implementation of [Flux Potential Analysis (FPA)](3_FPA), which quantitatively integrates continuous gene expression data with the model at pathway- and reaction-level. 

We also provide two auxiliary scriots. First is a python-based metabolic distance calculator that generates the shortest path from a given reaction to any other reaction in a given metabolic network. The [metabolic distance calculator](MetabolicDistance) provides the required input for applying FPA to any metabolic network. Second is CatExp, a collection of functions to categorize genes in a dataset into highly, moderatelt, lowly, and rarely expressed genes. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The main programs were developed and tested in MATLAB R2019a/R2017a. [COnstraint-Based Reconstruction and Analysis (The COBRA Toolbox)](https://opencobra.github.io/cobratoolbox/stable/) is required to perform the analysis. Check [here](https://opencobra.github.io/cobratoolbox/stable/installation.html) for the installation guidance of COBRA Toolbox. The programs were tested for COBRA Toolbox - 2020 version, but should be compatible with an earlier version. 

The Linear Program (LP) and Mixed-Integer Linear Problem (MILP) solver used in the study was [gurobi](http://gurobi.com) 8.10. The built-in solver interface of COBRA Toolbox was used, so that we expect our program to also work with other supported solver in COBRA Toolbox. Please [see here](https://opencobra.github.io/cobratoolbox/stable/installation.html#solver-installation) for furter information about solver availability. 

The auxiliary scripts in Python language were developed and tested in Python 2.7. 

### Installing

This package doesn't require any installation or compiling. Please see the following section for reproducing the tissue expression integration result in this study, or perform the analysis on a desired metabolic network model and gene expression dataset. 

## Running the tests

The repo includes five independent modules, [iMAT++ algorithm](1_iMAT++), [original iMAT algorithm](2_iMAT), [Flux Potential Analysis (FPA)](3_FPA), [metabolic distance calculator](MetabolicDistance), and [categorizer of gene expression data](CatExp). Please see the instruction within each module for running a test.

The followings are descriptions on each module (folder) listed.

[1_iMAT++](1_iMAT++): The iMAT++ module. This folder contains all the input, scripts and example codes for running the iMAT++ for <i>C. elegans</i> tissue-level integration and application on other models/datasets. 

[2_iMAT](2_iMAT): The original iMAT module. This folder contains all the input, scripts and example codes for running the original iMAT for <i>C. elegans</i> tissue-level integration. 

[3_FPA](3_FPA): The Flux Potential Analysis (FPA) module. This folder contains all the input, scripts and example codes for running the FPA for <i>C. elegans</i> tissue-level integration and application on other models/datasets. 

[MetabolicDistance](MetabolicDistance): The metabolic distance calculator module. This folder contains all the input, scripts and example codes for calculating the metabolic distance for any metabolic network model. The output distance matrix is required to run FPA. The output for the <i>C. elegans</i> models used in the reference study are available in the pertaining folders. 

[categorizer of gene expression data](CatExp): CatExp.py provides a set of functions for systematic categorization of gene expression levels. This folder also contains example inputs (gene expression datasets) and outputs (tables of categorized genes) used in the walkthrough help guide. The examples include the <i>C. elegans</i> tissue-level datasset used in the paper. 

[bins](bins): The shared functions for running the above mentioned analysis. These functions include modified version of some COBRA Toolbox functions and common new functions such as a molecular weight calculator.

[input](input): The shared input for running the above mentioned analysis. These inputs include <i>C. elegans</i> metabolic model and other input information.


## Contributing

Please contact us for reporting bugs or contributing purposes. Email contacting is encouraged. Please send to [Xuhang Li](mailto:xuhang.li@umassmed.edu) or [Safak Yilmaz](mailto:lutfu.yilmaz@umassmed.edu)


## Authors

* **Safak Yilmaz** - *Development of iMAT++, FPA, Metabolic Distance Calculator, and Gene Expression Categorizer* - [lsafak](https://github.com/lsafak)
* **Xuhang Li** - *Matlab implementation of iMAT++/iMAT/FPA and additional applications with a human model* - [XuhangLi](https://github.com/XuhangLi)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* TBD
