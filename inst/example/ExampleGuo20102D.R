# Example for running the tool on the Guo's Dataset
library(merlot)

# Read the example Guo dataset that is distributed with the package
DataFile= paste(find.package("merlot"), "/example/Guo2010.txt", sep="")
Dataset=ReadDataset(DataFile)

# Load the cell types
CellTypes=read.table(file=paste(find.package("merlot"), "/example/GuoFeatures.txt", sep=""), sep="\t", header = F, stringsAsFactors = F)
CellTypes=CellTypes[,2]
selected_colors=c("red", "orange", "yellow", "green", "cyan", "darkblue")

guo_colorcells=c()
guo_colorcells[which(CellTypes=="2C")]="red"
guo_colorcells[which(CellTypes=="4C")]="orange"
guo_colorcells[which(CellTypes=="8C")]="yellow"
guo_colorcells[which(CellTypes=="16C")]="green"
guo_colorcells[which(CellTypes=="32C")]="cyan"
guo_colorcells[which(CellTypes=="64C")]="darkblue"

# Embed Cells into their manifold
# library(destiny)
# DatasetDM <- DiffusionMap(Dataset$ExpressionMatrix, density.norm = T, verbose = F, sigma="global")
# CellCoordinates=DatasetDM@eigenvectors[,1:3]
# t <- -25
# theta = t / 180 * pi
# rot = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2, ncol=2)
# coords = DatasetDM@eigenvectors[,c(1,3)]
# tcoords = coords %*% rot
# # CellCoordinates=cbind(DatasetDM@eigenvectors[,2], tcoords[,1])

# The first 3 diffusion map components will be used for this example
CellCoordinates=read.table(file="/home/gonzalo/merlot/inst/example/GuoRotatedCoordinates.txt", sep="\t", header = F, stringsAsFactors = F)
CellCoordinates=as.matrix(CellCoordinates[,1:2])

# We calculate the scaffold tree using the first 3 diffusion components from the diffusion map
ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates, BranchMinLengthSensitive = sqrt(428))
# Plot the calculated tree
svg(filename = "/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Figuras/Fig1Bis/DMGuoScaffold.svg", height = 6, width = 6)
plot_scaffold_tree(ScaffoldTree = ScaffoldTree, colorcells = guo_colorcells)
dev.off()

NumberOfNodes=100
# We calculate the elastic principal tree using the scaffold tree for its initialization
ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes, FixEndpoints = F)

svg(filename = "/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Figuras/Fig1Bis/DMGuo.svg", height = 6, width = 6)
plot_elastic_tree(ElasticTree)
# plot_elastic_tree(ElasticTree, colorcells = guo_colorcells)
dev.off()

# Embedd the principal elastic tree on the gene expression space from which it was calculated.
EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = ElasticTree, increaseFactor_mu = 10, increaseFactor_lambda = 10)

# Calculate Pseudotimes for the nodes in the Tree in the full gene expression space
Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0=1)
plot_pseudotimes(CellCoordinates, Pseudotimes)

# Plot gene expression profile as a function of pseudotime
plot_pseudotime_expression_gene(GeneName = "Gata4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Runx1" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Nanog" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Fgf4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Fgfr2" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Klf2" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)

# Downstream analysis -------------
OrderedMatrix=plot_heatmaps_embedding(Pseudotimes, EmbeddedTree, log_tranform=F)

# Differentially Expressed Genes among two subpopulations in the tree
Group1=EmbeddedTree$Branches[[1]]
Group2=EmbeddedTree$Branches[[2]]
DifferentiallyExpressedGenes=subpopulations_differential_expression(SubPopulation1 = Group1, SubPopulation2 = Group2, EmbeddedTree = EmbeddedTree, mode = "cells")

# Differentially Expressed Genes in a specific branch
Branch1Genes=branch_differential_expression(Branch =1, EmbeddedTree, mode="cells")
Branch2Genes=branch_differential_expression(Branch =2, EmbeddedTree, mode="cells")

GetGeneCorrelationNetwork(EmbeddedTree$Nodes, cor_threshold = 0.7)
GetGeneCorrelationNetwork(Dataset$ExpressionMatrix, cor_threshold = 0.2)

# Differentially Expressed Genes among two subpopulations in the tree
Group1=EmbeddedTree$Branches[[4]]
Group2=EmbeddedTree$Branches[[5]]

DifferentiallyExpressedGenes=subpopulations_differential_expression(SubPopulation1 = Group1, SubPopulation2 = Group2, EmbeddedTree = EmbeddedTree, mode = "cells")

plot_pseudotime_expression_gene(GeneName = "Gata4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Runx1" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Fgf4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Fgfr2" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)


