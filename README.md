# MERLot: A MEthod for Reconstructing Lineage tree Topologies using scRNA-seq data

MERLot is a tool that can reconstruct the lineage tree topology that explains the emergence of different cell types from a progenitor population. MERLot is an R package than can reconstruct complex lineage tree topologies using coordinates for cells in a given manifold(like diffusion maps) as input.

## Get ready

## 1) R Dependencies:
MERLoT depends on certain R packages in order to work properly. Most of the packages can be installed either via CRAN (with the install.packages() function) or via Bioconductor.

* car
* rgl
* igraph
* fields
* ElPiGraph.r
* FactoMineR
* mclust

The Destiny package for creating diffusion maps was one of the dinmensionality reduction techniques we used in order to reconstruct lineage tree topologies in a low dimensional manifold.

The destiny package as well as instructions about how to install and use it can be found [here](https://www.helmholtz-muenchen.de/icb/research/groups/quantitative-single-cell-dynamics/software/destiny/index.html)

Optional packages:
* energy (needed for finding differentially expressed genes)
* VGAM

ElPiGraph.R can be installed following the instructions from the [developer's site](https://github.com/Albluca/ElPiGraph.R).

`install.packages("devtools")`
`library(devtools)`
`install.packages(c("bigpca", "irlba", "nsprcomp", "plotly","fields", "igraph", "rgl", "tictoc"))`


## 2) Installation

### a. Install from source
Either download the latest version from the git repository with `git clone https://github.com/soedinglab/merlot.git` or download one of our [releases](https://github.com/soedinglab/merlot/releases) and unpack it. Then, in R, call:

`install.packages("/path/to/merlot/folder",  types="source", repos = NULL)`

### b. Install directly from `git`

`library(devtools)`

`install_github("soedinglab/merlot")`

## 3) Known issues

#### Problems with `rgl`
The `rgl` library requires an X server to be present. Installing [xquartz](https://www.xquartz.org/) fixed this for us.