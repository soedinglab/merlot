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
library(destiny)
DatasetDM <- DiffusionMap(Dataset$ExpressionMatrix, density.norm = T, verbose = F, sigma="global")
#
t <- -20
theta = t / 180 * pi
rot = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2, ncol=2)
coords = DatasetDM@eigenvectors[,c(1,3)]
tcoords = coords %*% rot
plot(DatasetDM@eigenvectors[,2], tcoords[,1], main=t, pch=16, xlab="Component 2", ylab="Transformed Component", col=guo_colorcells)
# CellCoordinates=cbind(DatasetDM@eigenvectors[,2], tcoords[,1])
# write.table(CellCoordinates,file = "/home/gonzalo/Desktop/benchmark10/Guocoords.txt",  quote = F, row.names = F, col.names = F, sep="\t")

# The first 3 diffusion map components will be used for this example
CellCoordinates=DatasetDM@eigenvectors[,1:3]
# End Embedding into manifold



# We calculate the scaffold tree using the first 3 diffusion components from the diffusion map
ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates)
# Plot the calculated tree
plot_scaffold_tree(ScaffoldTree = ScaffoldTree, colorcells = guo_colorcells)

NumberOfNodes=100
# We calculate the elastic principal tree using the scaffold tree for its initialization
ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes, FixEndpoints = F)
plot_elastic_tree(ElasticTree, colorcells=guo_colorcells)

# Embedd the principal elastic tree on the gene expression space from which it was calculated.
EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = ElasticTree, increaseFactor_mu = 10, increaseFactor_lambda = 10)

# Calculate Pseudotimes for the nodes in the Tree in the full gene expression space
Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0=3)
plot_pseudotimes(CellCoordinates, Pseudotimes)

# Plot gene expression profile as a function of pseudotime
plot_pseudotime_expression_gene(GeneName = "Gata4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Fgf4" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Id2" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Sox2" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Nanog" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)
plot_pseudotime_expression_gene(GeneName = "Klf2" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)


# Rotation of coordinates for generating a 2D plot for this dataset
# rotate
t <- -20
theta = t / 180 * pi
rot = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2, ncol=2)
coords = DatasetDM@eigenvectors[,c(1,3)]
tcoords = coords %*% rot
plot(DatasetDM@eigenvectors[,2], tcoords[,1], main=t, pch=16, xlab="Component 2", ylab="Transformed Component")

CoordinatesToPlot=cbind(DatasetDM@eigenvectors[,2], tcoords[,1])
ElasticTreeToPlot=ElasticTree

coords = ElasticTreeToPlot$Nodes[,c(1,3)]
tcoords = coords %*% rot
plot(ElasticTreeToPlot$Nodes[,2], tcoords[,1], main=t, pch=16, xlab="Component 2", ylab="Transformed Component")
ElasticTreeToPlot$Nodes=cbind(ElasticTreeToPlot$Nodes[,2], tcoords[,1])
ElasticTreeToPlot$CellCoords=CoordinatesToPlot
plot_elastic_tree(ElasticTreeToPlot, colorcells=guo_colorcells)

# Plot the principal elastic tree
svg(filename = paste("/home/gonzalo/merlot/inst/example/Guo.svg", sep=""), height = 8, width = 8)
plot_elastic_tree(ElasticTree, legend = F)
dev.off()
# ElasticTreeInterp=DuplicateTreeNodes(ElasticTree)
# ElasticTreeInterp=DuplicateTreeNodes(ElasticTreeInterp)
# ElasticTreeInterp=DuplicateTreeNodes(ElasticTreeInterp)
# plot_elastic_tree(ElasticTreeInterp, colorcells = NULL)

# open3d()
# plot_elastic_tree(ElasticTree, colorcells = NULL)




svg(filename = "/home/gonzalo/Dropbox/SoedingGroup/ISMB2017/Show1.svg", height = 4, width = 6)
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

GetGeneCorrelationNetwork(EmbeddedTree$Nodes, cor_threshold = 0.7)
GetGeneCorrelationNetwork(Dataset$ExpressionMatrix, cor_threshold = 0.2)

#cluster into a K-Star with K=3
library(kbranches)
clust=kbranches.global(CellCoordinates, Kappa=3)

#plot cluster assignments
set.seed(1)
pal=RColorBrewer::brewer.pal(8,"Dark2")#change the color palette
plot(CellCoordinates,col=pal[clust$cluster],pch=21,bg=pal[clust$cluster])

# --------
# Load the cell types
CellTypes=read.table(file="/home/gonzalo/Desktop/Postdoc/Colabs/Niko/Papers/PreparingDatasets/ProcessedDatasets/Features/GuoFeatures.txt", sep="\t", header = F, stringsAsFactors = F)
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

#
# N_yk=80
# lambda_0=2.03e-09
# mu_0=0.00625
# # Calculate Tree Total Distance
# TreeTotalDistance=0
# dista_branch_i=c()
# for(i in 1:5)
# {
#   dista_branch_i=c(dista_branch_i, as.matrix(dist(CellCoordinates[ScaffoldTree$Branches[i,],], method = "euclidean", diag = FALSE, upper = TRUE, p = 2))[1,2])
#   TreeTotalDistance= TreeTotalDistance +dista_branch_i[i]
# }
#
# Branch_N_yk=c()
#
# plot3d(CellCoordinates)
# for(i in 1:5)
# {
#   ScaffoldTree$Branches[i,]
#
#   Branch_N_yk=c(Branch_N_yk, round(dista_branch_i[i]*N_yk/TreeTotalDistance))
#   # scaling according to N_yk
#   mu=(N_yk-1)*mu_0
#   lambda=((N_yk-2)**3)*lambda_0
#
#   ElasticTreePart <- computeElasticPrincipalGraph(Data = ScaffoldTree$CellCoordinates[ScaffoldTree$Cells2BranchesAssignments[[i]],], NumNodes = Branch_N_yk[i], NodesPositions = CellCoordinates[ScaffoldTree$Branches[i,],], Edges = c(1,2), Method = 'DefaultPrincipalTreeConfiguration', EP=lambda, RP=mu)
#   # ElasticTreePart[[1]]$Nodes[c(1,2),]=CellCoordinates[ScaffoldTree$Branches[i,],]
#   # points3d(CellCoordinates[ScaffoldTree$Branches[i,],], pch=16, col="green")
#   points3d(ElasticTreePart[[1]]$Nodes, pch=16, col=i, size=7)
#
# }
# plot3d(CellCoordinates)
# points3d(ElasticTree$Nodes, pch=16, col=i, size=7)
#
#
#
# ElasticTreePart <- computeElasticPrincipalGraph(Data = ScaffoldTree$CellCoordinates, NumNodes = 150, NodesPositions = ElasticTree$Nodes, Edges = ElasticTree$Edges, Method = 'CurveConfiguration', EP=lambda, RP=mu)
# plot3d(CellCoordinates)
# points3d(ElasticTreePart[[1]]$Nodes, pch=16, col=i, size=7)
