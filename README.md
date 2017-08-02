# MERLot: A MEthod for Reconstructing Lineage tree Topologies using scRNA-seq data

MERLot is a tool that can reconstruct the lineage tree topology that explains the emergence of different cell types from a progenitor population. MERLot is an R package than can reconstruct complex lineage tree topologies using coordinates for cells in a given manifold(like diffusion maps) as input.

## Get ready

## 1) Python Dependencies:
MERLoT consists of 1 part written in Python, which is distributed with the R package for which the following packages need to be installed. Take into account that MERLoT uses python 3.

* scipy
* pandas
* python3-tk
* numpy
* csgraph_mod (modified version of csgraph) which you can install from here: https://github.com/soedinglab/csgraph_mod

**NOTE:**
In case you don’t have a standard python3 installation, e.g you installed it using anaconda, when using the package you will need to set the location of your working python3 binary in the _python_url_ variable in the _ScaffoldTree()_ function. By default it is set to “/usr/bin/python3” (See the Vignette Section, ScaffoldTree() function, for more information).

## 2) R Dependencies:
MERLoT depends on certain R packages in order to work properly. Most of the packages can be installed either via CRAN (with the install.packages() function) or via Bioconductor.

* car
* rgl
* rpgraph
* igraph
* fields

The Destiny package for creating diffusion maps was one of the dinmensionality reduction techniques we used in order to reconstruct lineage tree topologies in a low dimensional manifold.

The destiny package as well as how to install it and use it can be found [here](https://www.helmholtz-muenchen.de/icb/research/groups/quantitative-single-cell-dynamics/software/destiny/index.html)

Optional packages:
* energy (needed for finding differentially expressed genes)
* VGAM

Rpgraph can be installed following the instructions from the [developer's site](https://github.com/Albluca/rpgraph/wiki).

The steps can be summarized in:
`install.packages(pkgs = "rJava", repos="http://rforge.net", type = 'source')`

For rJava **you have to have Java installed** in your system. You can install **default-jre, open jdk**

`install.packages("devtools")`
`library(devtools)`
`install.packages(c("bigpca", "irlba", "nsprcomp", "plotly","fields", "igraph", "rgl", "tictoc"))`


## 3) Download
Download the [merlot](https://github.com/soedinglab/merlot/archive/master.zip) file.


## 4) Install MERLoT
Install from source:

`install.packages("/url_in_your_computer/merlot-master.zip",  types="source", repos = NULL)`

Install from github:

`library(devtools)`

`install_github("soedinglab/merlot")`
