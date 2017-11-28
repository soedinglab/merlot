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

# The first 3 diffusion map components will be used for this example
CellCoordinates=DatasetDM@eigenvectors[,1:3]
# End Embedding into manifold

# We calculate the scaffold tree using the first 3 diffusion components from the diffusion map
ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates, BranchMinLengthSensitive = sqrt(428))
# Plot the calculated tree
plot_scaffold_tree(ScaffoldTree = ScaffoldTree, colorcells = guo_colorcells)
legend(x="bottomright", legend=c("2C", "4C", "8C", "16C", "32C", "64C"), col=selected_colors, pch=16, cex=0.7)

NumberOfNodes=100
# We calculate the elastic principal tree using the scaffold tree for its initialization
ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes)
plot_elastic_tree(ElasticTree, colorcells=guo_colorcells)

# Embedd the principal elastic tree on the gene expression space from which it was calculated.
EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = ElasticTree)

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

OrderedMatrix=plot_heatmaps_embedding(Pseudotimes, EmbeddedTree, log_tranform=F)

# Differentially Expressed Genes among two subpopulations in the tree
# Selecting cells in branch 1
Group1=EmbeddedTree$Branches[[1]]
# Selecting cells in branch 2
Group2=EmbeddedTree$Branches[[2]]

# Calculating whether each gene in the Expression Matrix is differentially expressed or not between the two selected populations
DifferentiallyExpressedGenes=subpopulations_differential_expression(SubPopulation1 = Group1, SubPopulation2 = Group2, EmbeddedTree = EmbeddedTree)

# Differentially Expressed Genes in a specific branch
Branch1Genes=branch_differential_expression(Branch =1, EmbeddedTree)
Branch2Genes=branch_differential_expression(Branch =2, EmbeddedTree)

# Calculating gene networks based on the original gene expression profiles
GetGeneCorrelationNetwork(Dataset$ExpressionMatrix, cor_threshold = 0.7)

# Calculating gene networks based on denoised expression profiles
GetGeneCorrelationNetwork(EmbeddedTree$Nodes, cor_threshold = 0.7)
