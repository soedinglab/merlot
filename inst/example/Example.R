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
EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = ElasticTree)

# Calculate Pseudotimes for the nodes in the Tree in the full gene expression space
Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0 = 1)

plot3d(CellCoordinates, col=col1)

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

mypal <- colorRampPalette( c( "yellow", "orange", "darkorange", "red", "black" ) )( length(Pseudotimes$Times_cells) )
col1=map2color(Pseudotimes$Times_cells, mypal)

# Requires WGCNA library to be loaded
plot_pseudotimes(CellCoordinates, Pseudotimes)

# Plot gene expression profile as a function of pseudotime
plot_pseudotime_expression_gene(GeneName = "Fgf4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)

png(filename = "/home/gonzalo/Desktop/merlot/WikiImages/GeneProfile.png", height = 400, width = 600)
plot_pseudotime_expression_gene(GeneName = "Gata4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
dev.off()
# Plots the averaged expression and the interpolated expression matrices. both matrices are returned as part of a list object
OrderedMatrix=plot_heatmaps_embedding(Pseudotimes, EmbeddedTree, log_tranform=F)

# Differentially Expressed Genes among two subpopulations in the tree
Group1=EmbeddedTree$Branches[[1]]
Group2=EmbeddedTree$Branches[[2]]
DifferentiallyExpressedGenes=subpopulations_differential_expression(SubPopulation1 = Group1, SubPopulation2 = Group2, EmbeddedTree = EmbeddedTree, mode = "cells")

# Differentially Expressed Genes in a specific branch
Branch1Genes=branch_differential_expression(Branch =1, EmbeddedTree, mode="cells")
Branch2Genes=branch_differential_expression(Branch =2, EmbeddedTree, mode="cells")
