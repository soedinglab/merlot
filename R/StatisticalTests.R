#' Subpopulations differential expression
#'
#' Given an elastic tree it calculates genes that are differentially expressed between two given populations in the tree. It applies a Kruskal-Wallis test.
#'
#' @param Subpopulation1 It is a vector containing the indexes of cells to be part of the first Subpopulation
#' @param Subpopulation2 It is a vector containing the indexes of cells to be part of the second Subpopulation
#' @param EmbeddedTree Embedded Elastic Tree as calculated by GenesSpaceEmbedding()
#' @param mode whether the cell expression profiles or the tree nodes expression profiles will be taken into account for the test. "cells" uses the expression values from the cells. "tree"
#'
#' @return genes_significances list object containing the p-values ($pvals) and e-values ($evals) vectors for the kruskall wallis test ordered from the highest significance to the lowest one. A vector containing the gene names ($GeneName) ordered as the significances vector is also provided. Another vector ($indexes) containing the original index order for genes in the columns from the expression matrix is also provided.
#'
#'
#' @export
subpopulations_differential_expression<-function(SubPopulation1, SubPopulation2,  EmbeddedTree, mode=c("tree", "cells"))
{
  # for testing
  # EmbeddedTree=EmbeddedTreeTopology
  # SubPopulation1= calc_path_elastic(2, 4, EmbeddedTreeTopology)
  # SubPopulation2= calc_path_elastic(3, 4, EmbeddedTreeTopology)
  # length(SubPopulation1)
  # length(SubPopulation2)
  # end for testing

  statistic_vector=c()

  # vectors for the expression values
  GenePopulation1=c()
  GenePopulation2=c()

  for (i in 1:dim(EmbeddedTree$Nodes)[2])
  {
    if(mode=="tree")
    {
      GenePopulation1=EmbeddedTree$Nodes[SubPopulation1, i]
      GenePopulation2=EmbeddedTree$Nodes[SubPopulation2, i]
    }
    else if(mode=="cells")
    {
      # take cells that are mapped to the nodes in populations1 and 2
      MappedPopulation1= EmbeddedTree$Cells2TreeNodes[which(EmbeddedTree$Cells2TreeNodes[,2] %in% SubPopulation1), 1]
      MappedPopulation2= EmbeddedTree$Cells2TreeNodes[which(EmbeddedTree$Cells2TreeNodes[,2] %in% SubPopulation2), 1]
      GenePopulation1=EmbeddedTree$CellCoords[MappedPopulation1, i]
      GenePopulation2=EmbeddedTree$CellCoords[MappedPopulation2, i]
    }

    if(sd(EmbeddedTree$Nodes[, i])==0)
    {
      # the gene has sd==0 and hence no differences can be observed
      statistic_vector=c(statistic_vector, 200)
    }
    else
    {
      BothPopulations <- data.frame(Y=c(GenePopulation1, GenePopulation2), group=c(rep("Population1", length(GenePopulation1)), rep("Population2", length(GenePopulation2))))
      KruskalResult=kruskal.test(Y~group, data=BothPopulations)
      if(!is.na(KruskalResult$p.value))
      {
        statistic_vector=c(statistic_vector, KruskalResult$p.value)
      }
      else
      {
        # the populations have expression vectors composed of zeros
        statistic_vector=c(statistic_vector, 300)
      }

    }
  }

  genes_significances=sort(statistic_vector, index.return=T)
  names(genes_significances)[names(genes_significances)=="x"] <- "pvals"
  names(genes_significances)[names(genes_significances)=="ix"] <- "indexes"
  genes_significances$GeneName <- colnames(EmbeddedTree$Nodes)[genes_significances$indexes]
  genes_significances$evals <- genes_significances$pvals * length(genes_significances$pvals)

  return(genes_significances)
}

#' differential expression in specific branch
#'
#' Given an elastic tree it calculates genes that are differentially expressed in a particular branch in the tree
#'
#' @param Branch It is a vector containing the indexes of cells to be part of the first Subpopulation
#' @param EmbeddedTree Embedded Elastic Tree as calculated by GenesSpaceEmbedding()
#' @param mode whether the cell expression profiles or the tree nodes expression profiles will be taken into account for the test. "cells" uses the expression values from the cells. "tree"
#'
#' @return genes_significances list object containing the p-values ($pvals) and e-values ($evals) vectors for the kruskall wallis test ordered from the highest significance to the lowest one. A vector containing the gene names ($GeneName) ordered as the significances vector is also provided. Another vector ($indexes) containing the original index order for genes in the columns from the expression matrix is also provided.
#'
#'
#' @export
branch_differential_expression<-function(Branch, EmbeddedTree, mode=c("tree", "cells"))
{
  BranchCells=EmbeddedTree$Branches[[Branch]]
  NumberNodes=dim(EmbeddedTree$Nodes)[1]
  NOT_Branch= which(!(seq(from=1, to=NumberNodes, by=1) %in% BranchCells))
  Branch_Genes=DifferentiallyExpressedGenes=subpopulations_differential_expression(SubPopulation1 = BranchCells, SubPopulation2 = NOT_Branch, EmbeddedTree = EmbeddedTree, mode = "cells")

  return(Branch_Genes)

}

#' Get Gene Correlation Network
#'
#' Given an gene expression matrix (raw or embedded in the high dimensional elastic tree) it calculates the correlations in between all pairs of gene expression profiles and plots the network
#'
#' @param ExpressionMatrix contaings the expression profiles(columns) along the cells (rows)
#' @param cor_threshold correlation threshold for showing edges in the network plot (default=0.4)
#' @param plot whether the network should be plotted or not (default=T)
#'
#' @return cor.genes pair correlation matrix
#'
#'
#' @export
GetGeneCorrelationNetwork <-function(ExpressionMatrix, cor_threshold=0.4, plot=T)
{
  GeneCorrelations<- cor(ExpressionMatrix, use="pair")
  GeneCorrelationsPlot=GeneCorrelations
  GeneCorrelationsPlot[GeneCorrelationsPlot < cor_threshold] <- 0
  diag(GeneCorrelationsPlot) <- 0
  graph <- graph.adjacency(GeneCorrelationsPlot, weighted=TRUE, mode="lower")
  set.seed(2)
  if(plot)
  {
    plot.igraph(graph, vertex.size=6, edge.width=0.5, vertex.label=rownames(GeneCorrelationsPlot), layout=layout.fruchterman.reingold(graph))
  }
  return(GeneCorrelations)
}



#' Filter Differentially Expressed Genes
#'
#' Given an expression matrix it keeps the differentially expressed genes in the data. First, genes expressed in less than 5% of the cells are filtered out.
#' Then a PCA is applied on the data. Next (optionally) a tSNE reduction is applied using the PCs components that were kept from the step before.
#' Clusters are foudn in the data, using the mclust package. The number of clusters is defined by using the BIC method.
#' Finally a kruskal.test is iteratively applied to find the differentially expressed genes in the data (e-val > 0.001)
#'
#' @param ExpressionMatrix contains the expression profiles(columns) along the cells (rows)
#' @param NPCs Number of principal components to be kept for clustering the data.
#'
#' @return DE_genes list of differentially expressed genes
#'
#'
#' @export
#'
FilterGenes<-function(ExprMatrix, NPCs=10)
{
  # testing
  # ExprMatrix=Dataset$ExpressionMatrix
  # end testing

  keep_genes=c()
  # Filter Genes without expression in at least 5% of cells
  for(i in 1:dim(ExprMatrix)[2])
  {
    if(length(which(ExprMatrix[,i]==0)) < dim(ExprMatrix)[1]*0.05 )
    {
      keep_genes=c(keep_genes, i)
    }
  }

  ExprMatrix=ExprMatrix[,keep_genes]
  # Calculate PCA for the data, keeping the first 50 componentes
  PCATransform=PCA(log2(t(ExprMatrix)+1), ncp = NPCs, graph = F)
  # Select those components that capture 60% of the variability
  # PCsToUse=which(PCATransform$eig[,3]<60)
  PCsToUse=1:NPCs
  ReducedSpace=PCATransform$var$coord[,PCsToUse]

  library(tsne)
  ReducedSpace=tsne(ReducedSpace)

  # library(Rtsne)
  # ReducedSpace=Rtsne(ReducedSpace, )
  # ReducedSpace=ReducedSpace$Y


  # We cluster the data with mclust
  library(mclust)
  # We decide the number of clusters to calculate with the BIC criteria
  mclusters=mclustBIC(ReducedSpace)

  # Cluster the data
  print("Clustering...")
  Clustering <- Mclust(ReducedSpace, x = mclusters)
  print(paste("Found ", length(unique(Clustering$classification)), "clusters.."))

  DE_genes=c()

  for(Cluster_i in 1:length(unique(Clustering$classification)))
  {
    print(paste("Finding DE genes in cluster", Cluster_i, sep=" "))
    DE_genes_clusteri=DE_genes_cluster(Cluster_i, Clustering, ExprMatrix)
    DE_genes=c(DE_genes, DE_genes_clusteri$GeneName[which(DE_genes_clusteri$evals<0.001)])
  }

  DE_genes=unique(DE_genes)

  return(DE_genes)
}

DE_genes_cluster <- function(Cluster, Clustering, ExprMatrix)
{
  # testing
  # Cluster=1
  # ExprMatrix=Dataset$ExpressionMatrix
  # end testing

  CellsClusteri=which(Clustering$classification==Cluster)
  AllOtherCells=which(!(1:dim(ExprMatrix)[1] %in% CellsClusteri))

  statistic_vector=c()

  for(Gene_j in 1:dim(ExprMatrix)[2])
  {
    # vectors for the expression values
    # Gene_j=1
    GenePopulation1=ExprMatrix[CellsClusteri, Gene_j]
    GenePopulation2=ExprMatrix[AllOtherCells, Gene_j]

    BothPopulations <- data.frame(Y=c(GenePopulation1, GenePopulation2), group=c(rep("Population1", length(GenePopulation1)), rep("Population2", length(GenePopulation2))))
    KruskalResult=kruskal.test(Y~group, data=BothPopulations)

    statistic_vector=c(statistic_vector, KruskalResult$p.value)
  }

  statistic_vector

  genes_significances=sort(statistic_vector, index.return=T)
  names(genes_significances)[names(genes_significances)=="x"] <- "pvals"
  names(genes_significances)[names(genes_significances)=="ix"] <- "indexes"
  genes_significances$GeneName <- colnames(ExprMatrix)[genes_significances$indexes]
  genes_significances$evals <- genes_significances$pvals * length(genes_significances$pvals)

  return(genes_significances)
}


