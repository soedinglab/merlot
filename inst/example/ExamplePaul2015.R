library(merlot)

JobName="Paul"
# Folder where results and temporal objects will be stored
JobFolder="/home/gonzalo/Desktop/Postdoc/Colabs/Niko/Papers/PreparingDatasets/TreeTopology_Parra2016/"
# Full location and name for the expression matrix file
DataFile="/home/gonzalo/Desktop/Postdoc/TreeTopology_Parra2016/ProcessedDatasets/PaulDescription.txt"


# colors for paul
Types=read.table(file = "/home/gonzalo/Desktop/Postdoc/Colabs/Niko/Papers/PreparingDatasets/ProcessedDatasets/Trials/Paul/MAP.csv", header=F, stringsAsFactors = F, sep=",")
paul_colorcells=c()

for(i in 1:19)
{
  paul_colorcells[which(Types[,2]==i)]=colors()[494+i]
}

CellCoordinates=read.table("/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Figuras/Fig1Bis/paul_ddr.csv", stringsAsFactors = F, header=T, row.names = 1, sep=",")
dim(CellCoordinates)
plot(CellCoordinates, pch=16, col=paul_colorcells)

Dataset=ReadDataset(DataFile)

ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates, NEndpoints = 3)
plot_scaffold_tree(ScaffoldTree = ScaffoldTree)

NumberOfNodes=70
# We calculate the elastic principal tree using the scaffold tree for its initialization
ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes, FixEndpoints = F, input="skeleton")
plot_elastic_tree(ElasticTree, colorcells = NULL)

ElasticTree <- computeElasticPrincipalGraph(Data = ScaffoldTree$CellCoordinates, NumNodes = 100,
                                            NodesPositions = ElasticTree$Nodes, Edges = ElasticTree$Edges,
                                            Method = 'CurveConfiguration', EP=0.00191062, RP=0.61875)

ElasticTree <- computeElasticPrincipalGraph(Data = ScaffoldTree$CellCoordinates, NumNodes = 200,
                                            NodesPositions = ElasticTree[[1]]$Nodes, Edges = ElasticTree[[1]]$Edges,
                                            Method = 'CurveConfiguration', EP=0.01575766, RP=1.24375)

plot(CellCoordinates, col="red")
points(ElasticTree[[1]]$Nodes, pch=16, col="black")

EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = ElasticTree, increaseFactor_mu = 10, increaseFactor_lambda=10)
Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0=3)
plot_pseudotimes(Coordinates, Pseudotimes)

plot_pseudotime_expression_gene("Klf1" , EmbeddedTree, Pseudotimes, addlegend = F, selectedcolors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue"))
plot_pseudotime_expression_gene("Cebpa" , EmbeddedTree, Pseudotimes, addlegend = F, selectedcolors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue"))
plot_pseudotime_expression_gene("Gata2" , EmbeddedTree, Pseudotimes, addlegend = F, selectedcolors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue"))

plot_gene_on_map("Klf1", CellCoordinates, Dataset$ExpressionMatrix)
plot_gene_on_map("Klf1", CellCoordinates, Pseudotimes$PseudoExpressionMatrix)
plot_gene_on_map("Cebpa", CellCoordinates, Pseudotimes$PseudoExpressionMatrix)


