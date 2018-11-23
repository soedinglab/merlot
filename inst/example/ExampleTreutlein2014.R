library(merlot)
# Read the example Guo dataset that is distributed with the package
DataFile= paste(find.package("merlot"), "/example/Treutlein2014.txt", sep="")
Dataset=ReadDataset(DataFile)

# Embed Cells into their manifold, in this case we use Diffusion Maps as calculated by Destiny
library(destiny)
DatasetDM <- DiffusionMap(Dataset$ExpressionMatrix)

CellCoordinates=DatasetDM@eigenvectors[,1:2]

ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates)

plot_scaffold_tree(ScaffoldTree = ScaffoldTree)

NumberOfNodes=100

# We calculate the elastic principal tree using the scaffold tree for its initialization
ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes)
plot_elastic_tree(ElasticTree)

EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = ElasticTree)

# Calculate Pseudotimes for the nodes in the Tree in the full gene expression space
Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0=1)
plot_pseudotimes(CellCoordinates, Pseudotimes)

plot_pseudotime_expression_gene(GeneName = "Marcks" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = T)

Dataset$Descriptions

# Load the cell types
CellInfo=read.table(file=paste(find.package("merlot"), "/example/TreutleinCellTypes.txt", sep=""), sep=" ", header=T, stringsAsFactors = F)
CellTypes=CellInfo[,1]
CellTimes=CellInfo[,2]

# selected colors for the cells
selected_colors=c("lightblue", "blue", "lightgreen", "skyblue", "seagreen", "darkblue", "cyan", "darkgreen")
# selected_colors=c("lightblue", "blue", "lightgreen", "violet", "magenta", "yellow", "green", "red")
treutlein_colorcells=c()
Types=sort(unique(CellTypes))
for(i in 1:length(Types))
{
  treutlein_colorcells[which(CellTypes==Types[i])]=selected_colors[i]
}

plot(CellCoordinates, col=treutlein_colorcells, pch=16)
legend(x="topright", legend=Types, col = selected_colors, pch=16)
box(lwd=2)
