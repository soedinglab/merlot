library(merlot)

# Read cell assignments provided by the dataset authors
CellTypes= paste(find.package("merlot"), "/example/PaulCellsMapping.csv", sep="")
Types=read.table(file =CellTypes, header=F, stringsAsFactors = F, sep=",")
paul_colorcells=c()
selected_colors=c("cyan1", "cyan4", "darkcyan",  "blue", "darkblue", "blue4", "navy", "darkgoldenrod", "darkgoldenrod1", "darkgoldenrod1", "gold", "bisque", "palegreen", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4", "darkolivegreen", "chartreuse4", "darkgreen")
for(i in 1:19)
{
  paul_colorcells[which(Types[,2]==i)]=selected_colors[i]
}

# Generate DDRTree coordinates using Monocle2
library(monocle)
ExpressionData=paste(find.package("merlot"), "/example/Paul2015.txt", sep="")
raw <- read.table(ExpressionData, sep="\t", header=T, row.names=1, stringsAsFactors = T)
exprs <- t(as.matrix(raw))
data <- newCellDataSet(exprs, expressionFamily = negbinomial.size(), lowerDetectionLimit = 1)

# run monocle ----
data <- estimateSizeFactors(data)
data <- estimateDispersions(data) # fails

data <- detectGenes(data, min_expr=0.1)

disp_table <- dispersionTable(data)
ordering_genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 2*dispersion_fit)$gene_id

data <- setOrderingFilter(data, ordering_genes)
data <- reduceDimension(data)
data <- orderCells(data, reverse=FALSE)

CellCoordinates <- t(reducedDimS(data))
plot(CellCoordinates)
# --------------------------------------------------

# We read the coordinates calculated by Monocle2
# CellCoordinates=read.table("/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Figuras/Fig1Bis/paul_ddr.csv", stringsAsFactors = F, header=T, row.names = 1, sep=",")
dim(CellCoordinates)
plot(CellCoordinates, pch=16, col=paul_colorcells)

ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates, NEndpoints = 3)

plot_scaffold_tree(ScaffoldTree = ScaffoldTree, colorcells = paul_colorcells)
legend(x="bottomright", legend=paste("Cluster ", 1:19 ), col=selected_colors, pch=16, cex=0.7)

NumberOfNodes=100
# We calculate the elastic principal tree using the scaffold tree for its initialization
#We start with a tree with 70 nodes (start_N_yk) and iteratively add 10 nodes (step_N_yk) until 100 nodes (N_yk) are reached
# ElasticTree=CalculateElasticTreeConstrained(ScaffoldTree=ScaffoldTree,  N_yk=NumberOfNodes, start_N_yk=70, step_N_yk=10)
ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes, ExtraScaffoldNodes = T)

plot_elastic_tree(ElasticTree, colorcells = paul_colorcells)
legend(x="bottomright", legend=paste("Cluster ", 1:19 ), col=selected_colors, pch=16, cex=0.7)

# We read the Expression Matrix for the Paul Dataset
Dataset=ReadDataset(ExpressionData)

# We embed the elastic tree on the expression matrix
EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = ElasticTree, increaseFactor_mu = 10, increaseFactor_lambda=10)

# we calculate the pseudotime
Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0=3)
plot_pseudotimes(CellCoordinates, Pseudotimes)

# imputed expression values as a function of pseudotime
plot_pseudotime_expression_gene("Klf1" , EmbeddedTree, Pseudotimes, addlegend = F, selectedcolors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue"))
plot_pseudotime_expression_gene("Cebpa" , EmbeddedTree, Pseudotimes, addlegend = F, selectedcolors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue"))
plot_pseudotime_expression_gene("Gata2" , EmbeddedTree, Pseudotimes, addlegend = F, selectedcolors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue"))


# Raw gene expression values mapped in the manifold space
plot_gene_on_map("Klf1", CellCoordinates, Dataset$ExpressionMatrix)
# Imputed gene expression values mapped in the manifold space
plot_gene_on_map("Klf1", CellCoordinates, Pseudotimes$PseudoExpressionMatrix)

