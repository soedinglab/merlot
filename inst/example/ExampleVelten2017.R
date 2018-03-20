library(merlot)

# ------Read coordinates and cell types----------------------
StemnetCoordinates=paste(find.package("merlot"), "/example/STEMNETCoords.txt", sep="")
Genes=read.table(file = StemnetCoordinates, stringsAsFactors = F, sep="\t")

CellCoordinates=Genes[,2:3]

# We read the cell types for the dataset and we prepare the vector with the cell colors
StemnetCellTypes=paste(find.package("merlot"), "/example/SteinmetzCellTypes.txt", sep="")
CellTypes=read.table(file = StemnetCellTypes, sep="\t", stringsAsFactors = F)
CellTypes=as.array(CellTypes[,1])

# We assigned different colors to the different cell types on each branch
branches_colors=c()
selected_colors=c("darkorchid", "yellow3", "magenta", "darkorange3", "forestgreen", "blue")
Types=sort(unique(CellTypes))
for(i in 1:length(Types))
{
  branches_colors[which(CellTypes==Types[i])]=selected_colors[i]
}

#------Calculate the Scaffold tree
ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates)

# svg(filename = "/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Figuras/Fig1Bis/DMSteinmetScaffold.svg", height = 6, width = 6)
plot_scaffold_tree(ScaffoldTree = ScaffoldTree, colorcells = branches_colors)
legend(x="bottomright", legend=c("B Cell", "Eo/Baso/Mast", "Ery", "Mk", "Mono/DC", "Neutro"), col=selected_colors, pch=16, cex=0.7)
# dev.off()

# Nicer values
lambda_0_i=0.80e-09
mu_0_i=0.00250

for(lambda_0_i in seq(0.30e-09, 0.90e-09, by=1e-10))
{
#----------------- Calculate the Elastic tree
NumberOfNodes=100
# ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes, FixEndpoints = F)
ElasticTree= CalculateElasticTreeNew(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes, NBranchScaffoldNodes=F, lambda_0=lambda_0_i, mu_0=mu_0_i)

# svg(filename = "/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Figuras/Fig1Bis/DMSteinmetzCellColors.svg", height = 6, width = 6)
plot_elastic_tree(ElasticTree, legend = F, colorcells = branches_colors)
legend(x="bottomright", legend=paste("Elasticity ", "lambda:", lambda_0_i, "mu:",mu_0_i  ), col=selected_colors, pch=16, cex=0.7)
# dev.off()
}

# Read the dataset to be embedded in the low dimensional manifold
ExpressionData=paste(find.package("merlot"), "/example/Velten2017.txt", sep="")
Dataset=ReadDataset(ExpressionData)
EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = Dataset$ExpressionMatrix, ElasticTree = ElasticTree)

# Calculate Pseudotime
Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0 =9)
plot_pseudotimes(CellCoordinates, Pseudotimes)

# Plot reconstructed gene expression profiles
# Names for genes are in the object Dataset$GeneNames
plot_pseudotime_expression_gene("KLF1 (ENSG00000105610)" , EmbeddedTree, Pseudotimes, addlegend = F, selectedcolors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue"))
plot_gene_on_map("KLF1 (ENSG00000105610)", CellCoordinates, Dataset$ExpressionMatrix)

