#' Plot Scaffold Tree
#'
#' Plots the scaffold tree on top of the distribution of cells in 2 or 3 coordinates.
#'
#'@param CellCoordinates coordinates for the cells. This are the same coordinates that were used in the Scaffold Tree Calculation
#'@export


plot_scaffold_tree <- function(ScaffoldTree, colorcells="gray", dims=dim(ScaffoldTree$CellCoordinates)[2])
{
  if(dims==2)
  {
    # Plot 2d Dijkstra's tree
    plot(ScaffoldTree$CellCoordinates[,1], ScaffoldTree$CellCoordinates[,2], pch=16, col=alpha(colorcells, 0.6), xlab="Component 1", ylab="Component 2", cex=1)
    for(i in 1:dim(ScaffoldTree$Branches)[1])
    {
      Path=calculate_path(ScaffoldTree$Branches[i,1], ScaffoldTree$Branches[i,2], ScaffoldTree$DijkstraPredecesors)
      lines(ScaffoldTree$CellCoordinates[Path,1], ScaffoldTree$CellCoordinates[Path,2], lwd=3, col="black")
      points(ScaffoldTree$CellCoordinates[ScaffoldTree$Endpoints,1], ScaffoldTree$CellCoordinates[ScaffoldTree$Endpoints,2], col="black", pch=16, cex=2)
      points(ScaffoldTree$CellCoordinates[ScaffoldTree$Branchpoints,1], ScaffoldTree$CellCoordinates[ScaffoldTree$Branchpoints,2], col="black", pch=16, cex=2)
      box(lwd=2)
    }
  }
  else if (dims==3)
  {
    plot3d(ScaffoldTree$CellCoordinates[ScaffoldTree$Endpoints,], col="green", size = 7)
    points3d(ScaffoldTree$CellCoordinates[ScaffoldTree$Branchpoints,], col="red", size = 7)
    legend3d(x="topright", legend=c("Endpoints", "Branchpoints"), col=c("green", "red"), pch=16)

    for(i in 1:dim(ScaffoldTree$Branches)[1])
    {
      lines3d(ScaffoldTree$CellCoordinates[calculate_path(ScaffoldTree$Branches[i,1], ScaffoldTree$Branches[i,2], ScaffoldTree$DijkstraPredecesors),], lwd=2)
    }

    # rgl.postscript("/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Dijkstra.svg","svg")
    open3d()
    plot3d(ScaffoldTree$CellCoordinates, type ="p", size=7, col = colorcells)
    for(i in 1:dim(ScaffoldTree$Branches)[1])
    {
      lines3d(ScaffoldTree$CellCoordinates[calculate_path(ScaffoldTree$Branches[i,1], ScaffoldTree$Branches[i,2], ScaffoldTree$DijkstraPredecesors),], lwd=2)
    }
  }
  else
  {
    print ("This function cannot work with more than 3 dimensions. See help")
  }
}

#' Plot Elastic Tree
#'
#' Plots an elastic tree on top of the cell coordinates on which it was calculated
#'
#' @param ElasticTree Principal Elastic Tree calculated with the CalculateElasticTreeFunction
#' @param colorcells vector of colors for the cells
#' @export

plot_elastic_tree <- function(ElasticTree, colorcells=NULL, legend=T)
{
  if(legend==T)
  {
    # change margins
    opar <- par(no.readonly=TRUE)
    par(xpd=T, mar=par()$mar+c(0,0,0,3))
  }

  # color cells by branch
  if(is.null(colorcells))
  {
    selected_colors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue", "aquamarine", "magenta", "brown", "gray")
    CellBranches=ElasticTree$Cells2Branches
    Types=sort(unique(CellBranches))

    for(i in 1:length(Types))
    {
      colorcells[which(CellBranches==Types[i])]=selected_colors[i]
    }
  }
  if(dim(ElasticTree$CellCoords)[2]==2)
  {
    # Plot 2d elastic Tree
    plot(ElasticTree$CellCoords[,1], ElasticTree$CellCoords[,2], pch=16, col=alpha(colorcells, 0.6), xlab="Component 1", ylab="Component 2", cex=1.5)
    for(i in 1:length(ElasticTree$Branches)[1])
    {
      points(ElasticTree$Nodes[1:length(ElasticTree$Topology$Endpoints),1], ElasticTree$Nodes[1:length(ElasticTree$Topology$Endpoints),2], col="black", pch=16, cex=2)
      points(ElasticTree$Nodes[(length(ElasticTree$Topology$Endpoints)+1):(length(ElasticTree$Topology$Endpoints)+length(ElasticTree$Topology$Branchpoints)),1], ElasticTree$Nodes[(length(ElasticTree$Topology$Endpoints)+1):(length(ElasticTree$Topology$Endpoints)+length(ElasticTree$Topology$Branchpoints)),2], col="black", pch=16, cex=2)
      lines(ElasticTree$Nodes[ElasticTree$Branches[[i]],1], ElasticTree$Nodes[ElasticTree$Branches[[i]],2], col="black", lwd=3)
      points(ElasticTree$Nodes[ElasticTree$Branches[[i]], 1], ElasticTree$Nodes[ElasticTree$Branches[[i]], 2], col="black", pch=16, cex=1)
      box(lwd=2)
    }
    if(legend==T)
    {
      legend("topright", inset=c(-0.12,0), legend=seq(1, length(ElasticTree$Branches), 1), pch=c(16), col = selected_colors[1:length(ElasticTree$Branches)], title="Branch")
      par(opar)
    }
  }
  # plot 3d elastic tree
  else if(dim(ElasticTree$CellCoords)[2]==3)
  {
    plot3d(ElasticTree$CellCoords, type ="p", size=7, col = colorcells)
    for(i in 1:length(ElasticTree$Branches))
    {
      points3d(ElasticTree$Nodes[ElasticTree$Branches[[i]],], size=7, lwd=3, col="black")
      lines3d(ElasticTree$Nodes[ElasticTree$Branches[[i]],], size=7, lwd=3, col="black")
    }
    legend3d("topright", legend=seq(1, length(ElasticTree$Branches), 1), pch=c(16), col = selected_colors[1:length(ElasticTree$Branches)], title="Branch")

    # rgl.postscript("/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Interpolated.svg","svg")
  }
  else
  {
    print ("This function cannot work with more than 3 dimensions. See help")
  }

}

#' Plot pseudotime expression gene
#'
#' Plot Gene Expression Profile as a function of pseudotime
#'
#'@param GeneName an element from the Datasets$GeneNames object.
#'@param EmbeddedTree High dimensional principal elastic tree.
#'@param Pseudotimes object calculated by CalculatePseudotimes()
#'
#'
#' @export

plot_pseudotime_expression_gene <- function (GeneName, EmbeddedTree, Pseudotimes, selectedcolors=rainbow(length(Pseudotimes$Branches)), branch_tags=c(), addlegend=F,  range_y=c())
{
  # GeneName: element from the GeneNames list in the Dataset object
  # EmbeddedTree: Embedded tree used to interpolate the genetic profiles for the tree nodes
  # Pseudotimes: Calculated for the tree nodes
  # selectedcolors: colors for the branches, by defaul the rainbow palette is used.

  # # Testing--------
  # GeneName="gene_60"
  selected_colors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue", "aquamarine", "magenta", "brown", "gray")

  # End testing
  # # Testing-------
  #
  Gene=which(colnames(EmbeddedTree$CellCoords)==GeneName)
  ExpressionGene=EmbeddedTree$Nodes[,Gene]
  ExpressionGeneOriginal=EmbeddedTree$CellCoords[,Gene]

  color_yk=Pseudotimes$Times_yk
  for(i in 1:length(Pseudotimes$Branches))
  {
    color_yk[Pseudotimes$Branches[[i]]]=selectedcolors[i]
  }

  # Mapping colors for cells and making them semi transparent
  CellColors=selectedcolors[Pseudotimes$Cells2Branches]
  CellColors <- rgb(t(col2rgb(CellColors)), alpha=50, maxColorValue=255)

  if(is.null(range_y))
  {
    range_y=range(ExpressionGene)
  }

  if(addlegend==T)
  {
    # change margins
    opar <- par(no.readonly=TRUE)
    par(xpd=F, mar=par()$mar+c(0,0,0,3))
  }

  # plot(Pseudotimes$Times_cells, ExpressionGeneOriginal, col=CellColors, pch=16, xlab="Pseudotime assignment", ylab="Expression Level", main=paste(GeneName, "Expression"), cex=1, ylim=range_y)
  plot(Pseudotimes$Times_cells, ExpressionGeneOriginal, col=CellColors, pch=16, xlab="Pseudotime", ylab="Expression Level", main=GeneName, cex=1, ylim=range_y)
  # points(Pseudotimes$Times_yk, ExpressionGene, pch=16, xlab="Pseudotime assignment", ylab="Expression Level", main=paste(GeneName, "Expression"), col=color_yk, cex=1.5)
  points(Pseudotimes$Times_yk[EmbeddedTree$Topology$Endpoints], ExpressionGene[EmbeddedTree$Topology$Endpoints], pch=16, col="black", cex=1.5)
  points(Pseudotimes$Times_yk[EmbeddedTree$Topology$Branchpoints], ExpressionGene[EmbeddedTree$Topology$Branchpoints], pch=16, col="black", cex=1.5)
  for(i in 1:length(Pseudotimes$Branches))
  {
    branch_i=Pseudotimes$Branches[[i]]
    aux=sort(Pseudotimes$Times_yk[branch_i], index.return=T)
    lines(aux$x, EmbeddedTree$Nodes[branch_i[aux$ix],Gene], col=selectedcolors[i], lwd=3)
  }
  # if(!is.null(branch_tags) & addlegend==T)
  # {
  #   legend(x="topright", legend = branch_tags, col=selectedcolors, pch=16)
  # }
  # else if (addlegend==T)
  # {
  #   legend(x="topright", legend = paste("branch", 1:length(Pseudotimes$Branches)), col=selectedcolors, pch=16)
  # }
  box(lwd=2)

  if(addlegend==T)
  {
    par(xpd=T)
    legend("topright", inset=c(-0.12,0), legend=seq(1, length(ElasticTree$Branches), 1), pch=c(16), col = selected_colors[1:length(ElasticTree$Branches)], title="Branch")
    par(opar)
  }


}

# Auxiliar function
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


#' Plot cells pseudotime
#'
#' Plot pseudotime for the cells in the given cell coordinatges
#'
#'@param CellCoordinates matrix containing up to 3 coordinates for cells in the manifold space
#'@param Pseudotimes object calculated by CalculatePseudotimes()
#'
#'
#' @export
plot_pseudotimes <- function (CellCoordinates, Pseudotimes)
{
  if(dim(CellCoordinates)[2]==3)
  {
    # paleta<-colorRampPalette(colors=c("yellow", "orange", "red", "black"))
    # color<-paleta(200)
    # col1<-numbers2colors(Pseudotimes$Times_cells , colors=color)
    mypal <- colorRampPalette( c( "yellow", "orange", "darkorange", "red", "black" ) )( length(Pseudotimes$Times_cells) )
    col1=map2color(Pseudotimes$Times_cells, mypal)
    plot3d(CellCoordinates[,1], CellCoordinates[,2], CellCoordinates[,3], col=col1, pch=16, size = 7, xlab = "DC1", ylab="DC2", zlab="DC3")
  }
  else if(dim(CellCoordinates)[2]==2)
  {
    paleta<-colorRampPalette(colors=c("yellow", "orange", "red", "black"))
    color<-paleta(200)
    col1<-numbers2colors(Pseudotimes$Times_cells , colors=color)
    plot(CellCoordinates[,1], CellCoordinates[,2], col=col1, pch=16, cex=1.5, ylab="Cooordinate 1", xlab="Coordinate 2")
  }
  else
  {
    print ("This function cannot plot pseudotime for datasets with more than 3 dimensions. If your dataset has more dimensions you may want to make a proyection from those into 2 or 3 dimensions in order to visualize their pseudotimes")
  }
}

#' Plot Expression Matrix Embedding
#'
#' Plots two heapmaps one for the averaged gene expression calculated by assigning cells to their closest tree nodes and one for the expression matrix after the elastic principal tree interpolation is performed. Both matrices are returned as part of a list object
#'
#' @param Pseudotimes pseudotimes object as calculated by CalculatePseudotimes()
#' @param EmbeddedTree elastic principal tree embedded in the full gene expression space as calculated by GenesSpaceEmbedding()
#' @param log_transform whether or not gene expression values should be log-transformed to improve visualization. By default it is equal to False.
#' @param cluster_genes wheter or not gene columns in the matrix will be clustered according to their expression values
#' @export
#'
plot_heatmaps_embedding <-function(Pseudotimes, EmbeddedTree, log_tranform=F, cluster_genes=T)
{

  # Testing
  # EmbeddedTree=EmbeddedTree
  # log_tranform=F
  # cluster_genes=T
  # end testing

  CellByBranchTimes=c()
  CheckPseudotime=c()

  for(i in 1:length(Pseudotimes$Branches))
  {
    branch_i=Pseudotimes$Branches[[i]]
    aux=sort(Pseudotimes$Times_yk[branch_i], index.return=T)
    CellByBranchTimes=c(CellByBranchTimes, branch_i[aux$ix])
  }

  CellByBranchTimes=unique(CellByBranchTimes)

  OrderedAverageNodes=EmbeddedTree$AveragedNodes[CellByBranchTimes,]
  OrderedInterpoaltedNodes=EmbeddedTree$Nodes[CellByBranchTimes,]

  if(log_tranform==T)
  {
    MaxNegValue=min(OrderedAverageNodes)
    # if(min(OrderedAverageNodes)<0)
    # {
    #   MaxNegValue=min(OrderedAverageNodes)
    # }
    OrderedAverageNodes=log(OrderedAverageNodes+MaxNegValue+1)

    MaxNegValue=min(OrderedInterpoaltedNodes)
    # if(min(OrderedInterpoaltedNodes)<0)
    # {
    #   MaxNegValue=min(OrderedAverageNodes)
    # }
    OrderedInterpoaltedNodes=log(OrderedInterpoaltedNodes+MaxNegValue+1)
  }

  if(cluster_genes==T)
  {
    # Cluster the matrixes by genes in the interpolated matrix
    OrderedInterpoaltedNodes.hc <- hclust(as.dist(1-cor(OrderedInterpoaltedNodes)))
    OrderedInterpoaltedNodes=OrderedInterpoaltedNodes[, OrderedInterpoaltedNodes.hc$labels[OrderedInterpoaltedNodes.hc$order]]
    OrderedAverageNodes=OrderedAverageNodes[, OrderedInterpoaltedNodes.hc$labels[OrderedInterpoaltedNodes.hc$order]]
  }


  my.colors = colorRampPalette(c("black", "red","orange", "yellow", "white"))
  range_heatmap= c(min(c(range(OrderedAverageNodes)[1], range(OrderedInterpoaltedNodes)[1])), max(c(range(OrderedAverageNodes)[2], range(OrderedInterpoaltedNodes)[2])))
  # svg(filename = paste("/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Figuras/Fig2/", "Yk_averaged", ".svg", sep=""), width = 6, height = 3.5)
  image.plot(t(OrderedAverageNodes), xlab="Genes", ylab="Cells", col =my.colors(300) , zlim=range_heatmap, xaxt='n', yaxt='n')
  # dev.off()
  # svg(filename = paste("/home/gonzalo/Dropbox/SoedingGroup/PaperTree/Figuras/Fig2/", "Yk_interpolated", ".svg", sep=""), width = 6, height = 3.5)
  image.plot(t(OrderedInterpoaltedNodes), xlab="Genes", ylab="Cells",  col =my.colors(300), zlim=range_heatmap, xaxt='n', yaxt='n')
  # dev.off()

  EmbeddedMatrices= list(AveragedMatrix=OrderedAverageNodes, InterpolatedMatrix=OrderedInterpoaltedNodes)


  return(EmbeddedMatrices)
}


#' Plots Gene Expression on Cells Coordinates
#'
#' Plots the expression of a given gene on a given set of coordinates
#'
#' @param GeneName One of the elements in the Dataset$GeneNames object.
#' @param CellCoordinates Matrix with 2 or 3 coordinates per cell.
#' @param Dataset Dataset as generated by the ReadDataset() function
#' @param Average wheter or not expression for cells will be averaged over the nearest knn neighbors. Default Average=FALSE
#' @param knn number of neighbors over which expression for the cells will be averaged in case Average=TRUE. Default: knn=5
#' @param Palette colors that will be used for the cells by default= Palette=colorRampPalette(colors=c("yellow", "orange", "red","black")))
#' @export
#'
plot_gene_on_map <-function(GeneName,
                            CellCoordinates,
                            Dataset,
                            CoordLabels=c("Component 1", "Component 2", "Component 3"),
                            Average=F,
                            knn=5,
                            Palette=colorRampPalette(colors=c("yellow", "orange", "red","black")))
{
  if (!dim(CellCoordinates)[2] %in% c(2,3)) {
    stop("Dimensions of CellCoordinates need to be 2 or 3.")
  }

  paleta = Palette
  color = paleta(400)
  Gene = which(colnames(Dataset$ExpressionMatrix)==GeneName)
  ExpressionGene = as.numeric(Dataset$ExpressionMatrix[,Gene])

  # In order to color the range properly, all values are applied an offset to be equal or greater to 0
  if(min(ExpressionGene)<0){ExpressionGene=ExpressionGene+min(ExpressionGene*-1)}

  if(Average==T)
  {
    ExpressionGeneAveraged=c()
    knn=5
    Distances=as.matrix(dist(CellCoordinates, method = "euclidean", diag = FALSE, upper = TRUE, p = 2))

    for (i in  1:length(ExpressionGene))
    {
      ordered=sort(Distances[i,], index.return=T)
      names(ordered)[names(ordered)=="x"] <- "pvals"
      names(ordered)[names(ordered)=="ix"] <- "indexes"
      ExpressionGeneAveraged= c(ExpressionGeneAveraged, mean(ExpressionGene[ordered$indexes[1:knn]]))
    }
    col1<-map2color(ExpressionGeneAveraged, color)
    if (dim(CellCoordinates)[2] == 2) {
      plot(CellCoordinates, main=paste("Expression Gene: ", GeneName), xlab=CoordLabels[1], ylab=CoordLabels[2], xlim=range(CellCoordinates[,1]), ylim=range(CellCoordinates[,2]), col=col1, pch=16, xaxt='n', yaxt='n')
    }
    else {
      plot3d(CellCoordinates, main=paste("Expression Gene: ", GeneName), xlab=CoordLabels[1], ylab=CoordLabels[2], zlab=CoordLabels[3], xlim=range(CellCoordinates[,1]), ylim=range(CellCoordinates[,2]), zlim=range(CellCoordinates[,3]), col=col1, pch=16, size=7, xaxt='n', yaxt='n')
    }
  }
  else
  {
    col1<-map2color(ExpressionGene, color)
    if (dim(CellCoordinates)[2] == 2) {
      plot(CellCoordinates, main=paste("Expression Gene: ", GeneName), xlab=CoordLabels[1], ylab=CoordLabels[2], xlim=range(CellCoordinates[,1]), ylim=range(CellCoordinates[,2]), col=col1, pch=16, xaxt='n', yaxt='n')
    }
    else {
      plot3d(CellCoordinates, main=paste("Expression Gene: ", GeneName), xlab=CoordLabels[1], ylab=CoordLabels[2], zlab=CoordLabels[3], xlim=range(CellCoordinates[,1]), ylim=range(CellCoordinates[,2]), zlim=range(CellCoordinates[,3]), col=col1, pch=16, size=7, xaxt='n', yaxt='n')
    }
  }
}
