# Predicting tissue-relevant metabolism at system-, pathway-, reaction- and metabolite-levels using single-cell RNA-seq data
# this tag indicates this repo is under development
We developped the package, MERGE, to systematically integrate gene expression data from single cell RNA-seq with the metabolic network reconstruction of C. elegans, iCEL1314. MERGE is applicable to other expression datasets and other metabolic network models. 

The package includes three modules: (1) matlab implementation of [iMAT++ algorithm](1_iMAT++), (2) matlab implementation of [original iMAT algorithm](2_iMAT) (for reproducing purpose), and (3) matlab implementation of [Flux Potential Analysis (FPA)](3_FPA). Additionally, we provide a python-based metabolic distance calculator to calculate the metabolic distance in any given metabolic network reconstruction. The [metabolic distance calculator](MetabolicDistance) will generate the required input for applying FPA to any metabolic network. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The programs were developed and tested in MATLAB R2019a/R2017a. [COnstraint-Based Reconstruction and Analysis (The COBRA Toolbox)](https://opencobra.github.io/cobratoolbox/stable/) is required to perform the analysis. Check [here](https://opencobra.github.io/cobratoolbox/stable/installation.html) for the installation guidance of COBRA Toolbox. The programs were tested for COBRA Toolbox - 2020 version, but should be compatible with an earlier version. 

The Linear Program (LP) and Mixed-Integer Linear Problem (MILP) solver used in the study was [gurobi](http://gurobi.com) 8.10. The built-in solver interface of COBRA Toolbox was used, so that we expect our program to also work with other supported solver in COBRA Toolbox. Please [see here](https://opencobra.github.io/cobratoolbox/stable/installation.html#solver-installation) for furter information about solver availability. 

### Installing

This package doesn't require any installation or compiling. Please see the following section for reproducing the tissue expression integration result in this study, or perform the analysis on a desired metabolic network reconstruction. 

## Running the tests

The repo includes four independent modules, [iMAT++ algorithm](1_iMAT++), [original iMAT algorithm](2_iMAT), [Flux Potential Analysis (FPA)](3_FPA), and [metabolic distance calculator](MetabolicDistance). Please see the instruction within each module for running a test.

The followings are descriptions on each module (folder) listed.

[1_iMAT++](1_iMAT++): The iMAT++ module. This folder contains all the input, scripts and example codes for running the iMAT++ for C. elegans Tissue Expression integration and application on other models/datasets. 

[2_iMAT](2_iMAT): The original iMAT module. This folder contains all the input, scripts and example codes for running the original iMAT for C. elegans Tissue Expression integration. 

[3_FPA](3_FPA): The Flux Potential Analysis (FPA) module. This folder contains all the input, scripts and example codes for running the FPA for C. elegans Tissue Expression integration and application on other models/datasets. 

[MetabolicDistance](MetabolicDistance): The metabolic distance calculator module. This folder contains all the input, scripts and example codes for calculating the metabolic distance for any metabolic network reconstruction. The output distance matrix is required to run FPA on a non-C. elegans model.

[bins](bins): The shared functions for running the above mentioned analysis. These functions include modified version of some COBRA Toolbox functions and common new functions such as molecular weight calculator.

[input](input): The shared input for running the above mentioned analysis. These inputs include C. elegans metabolic model and other input informations.


## Contributing

Please contact us for reporting bugs or contributing purposes. Email contacting is encouraged. Please send to [Xuhang Li](mailto:xuhang.li@umassmed.edu) or [Safak Yilmaz](mailto:lutfu.yilmaz@umassmed.edu)


## Authors

* **Safak Yilmaz** - *Development of iMAT++ and FPA and Metabolic Distance Calculator* - [lsafak](https://github.com/lsafak)
* **Xuhang Li** - *Matlab implementation of iMAT++/iMAT/FPA* - [XuhangLi](https://github.com/XuhangLi)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* TBD
