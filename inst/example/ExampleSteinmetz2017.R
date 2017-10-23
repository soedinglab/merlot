library(merlot)

# ------Read coordinates and cell types----------------------
Genes=read.table(file = "/home/gonzalo/Downloads/Steinmetz/STEMNETCoordinates", stringsAsFactors = F, sep="\t")

CellCoordinates=Genes[,2:3]

# We read the cell types for the dataset and we prepare the vector with the cell colors
CellTypes=read.table("/home/gonzalo/merlot/inst/example/SteinmetzCellTypes.txt", sep="\t", stringsAsFactors = F)
CellTypes=as.array(CellTypes[,1])
branches_colors=c()
# selected_colors=c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue")
selected_colors=c("darkorchid", "yellow3", "magenta", "darkorange3", "forestgreen", "blue")
Types=sort(unique(CellTypes))
for(i in 1:length(Types))
{
  branches_colors[which(CellTypes==Types[i])]=selected_colors[i]
}

#------Calculate the Scaffold tree
ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates)

svg(filename = "/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Figuras/Fig1Bis/DMSteinmetScaffold.svg", height = 6, width = 6)
plot_scaffold_tree(ScaffoldTree = ScaffoldTree, colorcells = branches_colors)
legend(x="bottomright", legend=c("B Cell", "Eo/Baso/Mast", "Ery", "Mk", "Mono/DC", "Neutro"), col=selected_colors, pch=16, cex=0.7)
dev.off()

#----------------- Calculate the Elastic tree
NumberOfNodes=100
ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes, FixEndpoints = F)

svg(filename = "/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Figuras/Fig1Bis/DMSteinmetzCellColors.svg", height = 6, width = 6)
plot_elastic_tree(ElasticTree, legend = F, colorcells = branches_colors)
dev.off()

# Read the dataset to be embedded in the low dimensional manifold
Dataset=ReadDataset("/home/gonzalo/merlot/inst/example/Steinmetz2017.txt")
EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = ElasticTree, increaseFactor_mu = 10, increaseFactor_lambda = 10)

# Calculate Pseudotime
Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0 =9)
plot_pseudotimes(CellCoordinates, Pseudotimes)

# Plot reconstructed gene expression profiles
plot_pseudotime_expression_gene("KLF1 (ENSG00000105610)" , EmbeddedTree, Pseudotimes, addlegend = F, selectedcolors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue"))
plot_gene_on_map("KLF1 (ENSG00000105610)", CellCoordinates, Dataset$ExpressionMatrix)
dim(Dataset$ExpressionMatrix)

# svg(filename = "/home/gonzalo/Desktop/JakobPoster/KLF1.svg", height = 4, width = 6)
# plot_pseudotime_expression_gene(GeneName = "KLF1 (ENSG00000105610)" , EmbeddedTree = EmbeddedTree, Pseudotimes = Pseudotimes, addlegend = F, range_y = "tree")
# dev.off()
