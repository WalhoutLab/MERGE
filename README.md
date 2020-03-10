# Systematic Integration of Gene Expression Data to Understand Tissue-Level Compartmentalization of Metabolic Network Function in C. elegans
# this tag indicates the README file and repo are under development
We developped the package to systematically integrate single cell (also applicable to bulk) RNA-seq data to the metabolic network reconstruction of C. elegans. The package includes three parts: (1) matlab implantation of [iMAT++ algorithm](1_iMAT++), (2) matlab implatation of [original iMAT algorithm](2_iMAT) (for reproducing purpose), and (3) matlab implantation of [Flux Potential Analysis (FPA)](3_FPA). Additionally, we provide a python-based metabolic distance calculator to calculate the metabolic distance in any given metabolic network reconstruction. The [metabolic distance calculator](MetabolicDistance) will generate the required input for applying FPA to any metabolic network. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The programs were developed and tested in MATLAB R2019a/R2017a. [COnstraint-Based Reconstruction and Analysis (The COBRA Toolbox)](https://opencobra.github.io/cobratoolbox/stable/) is required to perform the analysis. Check [here](https://opencobra.github.io/cobratoolbox/stable/installation.html) for the installation guidance of COBRA Toolbox. The programs were tested for COBRA Toolbox - 2020 version, but should be compatible with an earlier version. 

The Linear Program (LP) and Mixed-Integer Linear Problem (MILP) solver used in the study was [gurobi]() 8.10. The built-in solver interface of COBRA Toolbox was used, so that we expect our program to also work with other supported solver in COBRA Toolbox. Please [see here](https://opencobra.github.io/cobratoolbox/stable/installation.html#solver-installation) for furter information about solver availability. 

### Installing

This package doesn't require any installation or compiling. Please see the following section for reproducing the tissue expression integration result in this study, or perform the analysis on a desired metabolic network reconstruction. 

## Running the tests

The repo includes four independent modules, [iMAT++ algorithm](1_iMAT++), [original iMAT algorithm](2_iMAT), [Flux Potential Analysis (FPA)](3_FPA), and [metabolic distance calculator](MetabolicDistance). Please see the instruction within each module for running a test.

The followings are descriptions on each folder listed.
[1_iMAT++](1_iMAT++): The iMAT++ module. This folder contains all the input, scripts and example codes for running the iMAT++ for C. elegans Tissue Expression integration. 
[2_iMAT](2_iMAT): The original iMAT module. This folder contains all the input, scripts and example codes for running the original iMAT for C. elegans Tissue Expression integration. 
[3_FPA](3_FPA): The Flux Potential Analysis (FPA) module. This folder contains all the input, scripts and example codes for running the FPA for C. elegans Tissue Expression integration. 
[MetabolicDistance](MetabolicDistance): The metabolic distance calculator module. This folder contains all the input, scripts and example codes for calculating the metabolic distance for any metabolic network reconstruction. The output distance matrix is required to run FPA on a non-C. elegans model.
[bins](bins): The shared functions for running the above mentioned analysis. These functions include modified version of a COBRA Toolbox function and common new functions such as molecular weight calculator.
[input](input): The shared input for running the above mentioned analysis. These inputs include C. elegans metabolic model and other input informations.


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.


## Authors

* **Safak Yilmaz** - *Development of iMAT++ and FPA and Metabolic Distance Calculator* - [PurpleBooth](https://github.com/lsafak)
* **Xuhang Li** - *Matlab Implantation of iMAT++/iMAT/FPA* - [PurpleBooth](https://github.com/XuhangLi)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* TBD
