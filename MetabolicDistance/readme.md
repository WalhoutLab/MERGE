This guidance shows how to use MetabolicDistance.py tool to find distances between reactions in a metabolic network model. 

### Reproducing the FPA results for the <i>C. elegans</i> dual tissue model

The following guidance, which can be run as a python program by copying the text in between the long lines, shows how to find distances from a particular reaction in the network to all other (reachable) reactions. See further below for other possible analyses with the metabolic distance tool, as well as options to use this tool efficiently in a high performance computing facility.



For other potential analyses, use help in other functions:
e.g.
>help(Net.findFrowardLoops)

To obtain a global distance matrix that shows the distance from every reaction (rows) to every other reaction (columns) the network must be traversed from every reaction using the above code. Doing this by placing the code in a for loop is an option, but can take many hours if not parallelized. How all distances can be efficiently calculated using a computer cluster is demonstrated using an example pipeline of three codes (see files with names formatted as "exampleCluster_xxx.py"). The example calculation was carried out in Massachusetts Green High Performance Computing Center (https://www.mghpcc.org/). This cluster uses LSF as job scheduler. The same code should be applicable in any high performance computing center using LSF, with minor modifications in job description if necessary. 



