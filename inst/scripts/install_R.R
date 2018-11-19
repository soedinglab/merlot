library("devtools")

# install ElPiGraph for the elastic principal trees
devtools::install_github("Albluca/distutils") # commit cd82bcac0aa919261813a542f2e4a4939eeecace
devtools::install_github("Albluca/ElPiGraph.R") # commit 251415810d1cf2a973eaa51d0f8e520e14dbe1cb

# isntall destiny for the diffusion maps
# (default dimensionality reduction)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("destiny", version = "3.8")

install.packages("fields")
install.packages(c("impute", "preprocessCore", "GO.db", "AnnotationDbi"))
install.packages("WGCNA")
install.packages("rgl")
install.packages("optparse")

devtools::install_github("soedinglab/merlot")