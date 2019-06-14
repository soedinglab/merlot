# Example for running the tool on the Guo's Dataset
library(merlot)

# Read the example Guo dataset that is distributed with the package
DataFile= paste(find.package("merlot"), "/example/GuoDescription.txt", sep="")
Dataset=ReadDataset(DataFile)

# Embed Cells into their manifold
library(destiny)
DatasetDM <- DiffusionMap(Dataset$ExpressionMatrix, density.norm = T, verbose = F, sigma="global")

# The first 3 diffusion map components will be used for this example
CellCoordinates=DatasetDM@eigenvectors[,1:3]
# End Embedding into manifold

# We calculate the scaffold tree using the first 3 diffusion components from the diffusion map
ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates)
# Plot the calculated tree
plot_scaffold_tree(ScaffoldTree = ScaffoldTree)

NumberOfNodes=100
# We calculate the elastic principal tree using the scaffold tree for its initialization
ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes, FixEndpoints = F)
# Plot the principal elastic tree
plot_elastic_tree(ElasticTree, colorcells = NULL)

# Embedd the principal elastic tree on the gene expression space from which it was calculated.
EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = ElasticTree, increaseFactor = 30)

# Calculate Pseudotimes for the nodes in the Tree in the full gene expression space
Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0=3)

# Requires WGCNA library to be loaded
plot_pseudotimes(CellCoordinates, Pseudotimes)

# Plot gene expression profile as a function of pseudotime
plot_pseudotime_expression_gene(GeneName = "Gata4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Fgf4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Id2" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Sox2" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Nanog" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Klf2" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)

# Plots the averaged expression and the interpolated expression matrices. both matrices are returned as part of a list object
OrderedMatrix=plot_heatmaps_embedding(Pseudotimes, EmbeddedTree, log_tranform=F)

# Differentially Expressed Genes among two subpopulations in the tree
Group1=EmbeddedTree$Branches[[1]]
Group2=EmbeddedTree$Branches[[2]]
DifferentiallyExpressedGenes=subpopulations_differential_expression(SubPopulation1 = Group1, SubPopulation2 = Group2, EmbeddedTree = EmbeddedTree, mode = "cells")

# Differentially Expressed Genes in a specific branch
Branch1Genes=branch_differential_expression(Branch =1, EmbeddedTree, mode="cells")
Branch2Genes=branch_differential_expression(Branch =2, EmbeddedTree, mode="cells")

#cluster into a K-Star with K=3
clust=kbranches.global(CellCoordinates,Kappa=3)

#plot cluster assignments
set.seed(1)
pal=RColorBrewer::brewer.pal(8,"Dark2")#change the color palette

plot(CellCoordinates,col=pal[clust$cluster],pch=21,bg=pal[clust$cluster])

# --------
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


DatasetDM <- DiffusionMap(Dataset$ExpressionMatrix, density.norm = T, verbose = F, sigma="global")

t <- -15
theta = t / 180 * pi
rot = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2, ncol=2)
coords = DatasetDM@eigenvectors[,c(1,3)]
tcoords = coords %*% rot
plot(DatasetDM@eigenvectors[,2], tcoords[,1], main=t, pch=16, xlab="Component 2", ylab="Transformed Component", col=guo_colorcells)
legend(x="bottomright", legend=c("2C", "4C", "8C", "16C", "32C", "64C"), col=selected_colors, pch=16)
CellCoordinatesRotated=cbind(DatasetDM@eigenvectors[,2], tcoords[,1])

# We calculate the scaffold tree using the first 3 diffusion components from the diffusion map
ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinatesRotated, NEndpoints = 4)
# Plot the calculated tree
plot_scaffold_tree(ScaffoldTree = ScaffoldTree)

NumberOfNodes=100
# We calculate the elastic principal tree using the scaffold tree for its initialization
ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes, FixEndpoints = F)
# Plot the principal elastic tree
plot_elastic_tree(ElasticTree, colorcells = NULL)

points(ElasticTree$Nodes[2,1], ElasticTree$Nodes[2,2], col="red", pch=16)

Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0=3)

plot_pseudotime_expression_gene(GeneName = "Gata4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Nanog" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Fgf4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Gata4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
