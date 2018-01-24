library(merlot)
# Read the example Guo dataset that is distributed with the package
DataFile= paste(find.package("merlot"), "/example/Treutlein2014.txt", sep="")
Dataset=ReadDataset(DataFile)

# Embed Cells into their manifold, in this case we use Diffusion Maps as calculated by Destiny
library(destiny)
DatasetDM <- DiffusionMap(Dataset$ExpressionMatrix, density.norm = T, verbose = F, sigma="global")

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
