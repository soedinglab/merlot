#' Plot Scaffold Tree
#'
#' Plots the scaffold tree on top of the distribution of cells in 2 or 3 coordinates.
#'
#' @param ScaffoldTree scaffoldTree calculated by the CalculateScaffoldTree function.
#' @param colorcells The colours to use for the cells
#' @param dims The number of dimensions to plot
#'
#' @export
#' @importFrom graphics plot box lines points legend
#' @importFrom rgl plot3d points3d legend3d lines3d open3d
#' @importFrom scales alpha
plot_scaffold_tree <- function(ScaffoldTree, colorcells="gray", dims=dim(ScaffoldTree$CellCoordinates)[2])
{
  if(dims==2)
  {
    # Plot 2d Dijkstra's tree
    plot(ScaffoldTree$CellCoordinates[,1], ScaffoldTree$CellCoordinates[,2], pch=16, col=scales::alpha(colorcells, 0.6), xlab="Component 1", ylab="Component 2", cex=1.5)

    if(length(ScaffoldTree$Endpoints)==2)
    {
      Path=calculate_path(ScaffoldTree$Branches[1], ScaffoldTree$Branches[2], ScaffoldTree$DijkstraPredecesors)
      lines(ScaffoldTree$CellCoordinates[Path,1], ScaffoldTree$CellCoordinates[Path,2], lwd=3, col="black")
      points(ScaffoldTree$CellCoordinates[ScaffoldTree$Endpoints,1], ScaffoldTree$CellCoordinates[ScaffoldTree$Endpoints,2], col="black", pch=16, cex=2)
      points(ScaffoldTree$CellCoordinates[ScaffoldTree$Branchpoints,1], ScaffoldTree$CellCoordinates[ScaffoldTree$Branchpoints,2], col="black", pch=16, cex=2)
      box(lwd=2)
    }else
    {
      for(i in 1:dim(ScaffoldTree$Branches)[1])
      {
        Path=calculate_path(ScaffoldTree$Branches[i,1], ScaffoldTree$Branches[i,2], ScaffoldTree$DijkstraPredecesors)
        lines(ScaffoldTree$CellCoordinates[Path,1], ScaffoldTree$CellCoordinates[Path,2], lwd=3, col="black")
        points(ScaffoldTree$CellCoordinates[ScaffoldTree$Endpoints,1], ScaffoldTree$CellCoordinates[ScaffoldTree$Endpoints,2], col="black", pch=16, cex=2)
        points(ScaffoldTree$CellCoordinates[ScaffoldTree$Branchpoints,1], ScaffoldTree$CellCoordinates[ScaffoldTree$Branchpoints,2], col="black", pch=16, cex=2)
        box(lwd=2)
      }
    }
  }
  else if (dims==3)
  {
    plot3d(ScaffoldTree$CellCoordinates[ScaffoldTree$Endpoints,], col="green", size = 7)
    points3d(ScaffoldTree$CellCoordinates[ScaffoldTree$Branchpoints,], col="red", size = 7)
    legend3d(x="topright", legend=c("Endpoints", "Branchpoints"), col=c("green", "red"), pch=16)

    if(length(ScaffoldTree$Endpoints)==2)
    {
      lines3d(ScaffoldTree$CellCoordinates[calculate_path(ScaffoldTree$Branches[1], ScaffoldTree$Branches[2], ScaffoldTree$DijkstraPredecesors),], lwd=2)
    }else
    {
      for(i in 1:dim(ScaffoldTree$Branches)[1])
      {
        lines3d(ScaffoldTree$CellCoordinates[calculate_path(ScaffoldTree$Branches[i,1], ScaffoldTree$Branches[i,2], ScaffoldTree$DijkstraPredecesors),], lwd=2)
      }
    }

    open3d()
    plot3d(ScaffoldTree$CellCoordinates, type ="p", size=7, col = colorcells)

    if(length(ScaffoldTree$Endpoints)==2)
    {
      lines3d(ScaffoldTree$CellCoordinates[calculate_path(ScaffoldTree$Branches[1], ScaffoldTree$Branches[2], ScaffoldTree$DijkstraPredecesors),], lwd=2)
    }else
    {
      for(i in 1:dim(ScaffoldTree$Branches)[1])
      {
        lines3d(ScaffoldTree$CellCoordinates[calculate_path(ScaffoldTree$Branches[i, 1], ScaffoldTree$Branches[i, 2], ScaffoldTree$DijkstraPredecesors),], lwd=2)
      }
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
#' @param legend Whether or not to plot the legend
#' @param legend_names The names to use for the legend
#' @export
#'
#' @importFrom graphics par plot points lines box legend
#' @importFrom rgl plot3d points3d legend3d lines3d open3d
#' @importFrom scales alpha
plot_elastic_tree <- function(ElasticTree, colorcells=NULL, legend=F, legend_names=c())
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
    selected_colors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue", "aquamarine", "magenta", "brown", "gray", "wheat1", "azure4", "lightsalmon4", "navy", "sienna1", "gold4", "red4", "violetred")
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
    plot(ElasticTree$CellCoords[,1], ElasticTree$CellCoords[,2], pch=16, col=scales::alpha(colorcells, 0.6), xlab="Component 1", ylab="Component 2", cex=1.5)
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
    plot3d(ElasticTree$CellCoords, type ="p", size=7, col = colorcells, xlab="Component 1", ylab="Component 2", zlab="Component 3")
    for(i in 1:length(ElasticTree$Branches))
    {
      points3d(ElasticTree$Nodes[ElasticTree$Branches[[i]],], size=7, lwd=3, col="black")
      lines3d(ElasticTree$Nodes[ElasticTree$Branches[[i]],], size=7, lwd=3, col="black")
    }

    if(legend==T)
    {
      if(is.null(legend_names))
      {
        legend3d("topright", legend=seq(1, length(ElasticTree$Branches), 1), pch=c(16), col = selected_colors[1:length(ElasticTree$Branches)], title="Branch")
      }
      else
      {
        legend3d("topright", legend=legend_names, pch=c(16), col = unique(colorcells), title="Branch")
      }

    }
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
#' @param GeneName an element from the Datasets$GeneNames object.
#' @param EmbeddedTree High dimensional principal elastic tree.
#' @param Pseudotimes object calculated by CalculatePseudotimes()
#' @param selectedcolors The colors to use for the branches
#' @param branch_tags The tags to use for the branches
#' @param addlegend Whether or not to add a legend
#' @param range_y Must be either "cells" or "tree"
#'
#'
#' @export
#' @importFrom grDevices rgb col2rgb rainbow
#' @importFrom graphics par plot points lines box legend
plot_pseudotime_expression_gene <- function (GeneName, EmbeddedTree, Pseudotimes, selectedcolors=grDevices::rainbow(length(Pseudotimes$Branches)), branch_tags=c(), addlegend=F,  range_y="tree")
{
  # GeneName: element from the GeneNames list in the Dataset object
  # EmbeddedTree: Embedded tree used to interpolate the genetic profiles for the tree nodes
  # Pseudotimes: Calculated for the tree nodes
  # selectedcolors: colors for the branches, by default the rainbow palette is used.

  # # Testing--------
  # GeneName="Fgf4"
  selected_colors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue", "aquamarine", "magenta", "brown", "gray", "wheat1", "azure4", "lightsalmon4", "navy", "sienna1", "gold4", "red4", "violetred")

  # End testing
  # # Testing-------
  #
  Gene=which(colnames(EmbeddedTree$CellCoords)==GeneName)
  ExpressionGene=EmbeddedTree$Nodes[,Gene]
  ExpressionGeneOriginal=EmbeddedTree$CellCoords[,Gene]

  color_yk=Pseudotimes$Times_yk
  for(i in 1:length(Pseudotimes$Branches))
  {
    color_yk[Pseudotimes$Branches[[i]]]=selected_colors[i]
  }

  # Mapping colors for cells and making them semi transparent
  CellColors=selected_colors[Pseudotimes$Cells2Branches]
  CellColors <- grDevices::rgb(t(grDevices::col2rgb(CellColors)), alpha=50, maxColorValue=255)

  if(range_y=="cells")
  {
    range_y=range(ExpressionGeneOriginal)
  }
  else if(range_y=="tree")
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

  if(length(EmbeddedTree$Topology$Endpoints)>2)
  {
    points(Pseudotimes$Times_yk[EmbeddedTree$Topology$Branchpoints], ExpressionGene[EmbeddedTree$Topology$Branchpoints], pch=16, col="black", cex=1.5)
  }

  for(i in 1:length(Pseudotimes$Branches))
  {
    # If and endpoint or branchpoint was selected as t0
    if(is.null(Pseudotimes$C0))
    {
      branch_i=Pseudotimes$Branches[[i]]
      aux=sort(Pseudotimes$Times_yk[branch_i], index.return=T)
      lines(aux$x, EmbeddedTree$Nodes[branch_i[aux$ix],Gene], col=selected_colors[i], lwd=3)
    }else{
      #If a cell was selected as t0
      if(Pseudotimes$T0 %in% Pseudotimes$Branches[[i]])
      {
        # if the branch contains t0, it will be plotted in two parts
        middle_node=which(Pseudotimes$Branches[[i]]==Pseudotimes$T0)
        semi_branch_i1=Pseudotimes$Branches[[i]][1:middle_node]
        aux=sort(Pseudotimes$Times_yk[semi_branch_i1], index.return=T)
        lines(aux$x, EmbeddedTree$Nodes[semi_branch_i1[aux$ix],Gene], col=selected_colors[i], lwd=3)

        semi_branch_i2=Pseudotimes$Branches[[i]][middle_node:length(Pseudotimes$Branches[[i]])]
        aux=sort(Pseudotimes$Times_yk[semi_branch_i2], index.return=T)
        lines(aux$x, EmbeddedTree$Nodes[semi_branch_i2[aux$ix],Gene], col=selected_colors[i], lwd=3)
      }else{
        # if branch does not contain t0 it is plotted as before
        branch_i=Pseudotimes$Branches[[i]]
        aux=sort(Pseudotimes$Times_yk[branch_i], index.return=T)
        lines(aux$x, EmbeddedTree$Nodes[branch_i[aux$ix],Gene], col=selected_colors[i], lwd=3)
      }
    }
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
    if(is.null(branch_tags))
    {
      legend("topright", inset=c(-0.12,0), legend=seq(1, length(EmbeddedTree$Branches), 1), pch=c(16), col = selected_colors[1:length(EmbeddedTree$Branches)], title="Branch")
    }
    else
    {
      legend("topright", inset=c(-0.16,0), legend=branch_tags, col = selected_colors[1:length(EmbeddedTree$Branches)], title="Branch", pch=c(16))
    }

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
#' @param CellCoordinates matrix containing up to 3 coordinates for cells in the manifold space
#' @param Pseudotimes object calculated by CalculatePseudotimes()
#' @param Projected boolean parameter to plot discrete suppport nodes pseudotime or projected pseudotime.
#'
#'
#' @export
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot legend
#' @importFrom rgl plot3d legend3d
plot_pseudotimes <- function (CellCoordinates, Pseudotimes, Projected=F)
{

  mypal <- grDevices::colorRampPalette( c( "yellow", "orange", "darkorange", "red", "black" ) )( length(Pseudotimes$Times_cells) )
  if(Projected==T)
  {
    col1=map2color(Pseudotimes$Projected_Times_Cells, mypal)
  }else
  {
    col1=map2color(Pseudotimes$Times_cells, mypal)
  }


  if(dim(CellCoordinates)[2]==3)
  {
    # paleta<-grDevices::colorRampPalette(colors=c("yellow", "orange", "red", "black"))
    # color<-paleta(200)
    # col1<-numbers2colors(Pseudotimes$Times_cells , colors=color)
    # mypal <- grDevices::colorRampPalette( c( "yellow", "orange", "darkorange", "red", "black" ) )( length(Pseudotimes$Times_cells) )
    # col1=map2color(Pseudotimes$Times_cells, mypal)
    plot3d(CellCoordinates[,1], CellCoordinates[,2], CellCoordinates[,3], col=col1, pch=16, size = 7, xlab = "DC1", ylab="DC2", zlab="DC3")
    legend3d()

  }
  else if(dim(CellCoordinates)[2]==2)
  {
    # mypal <- grDevices::colorRampPalette( c( "yellow", "orange", "darkorange", "red", "black" ) )( length(Pseudotimes$Times_cells) )
    # col1=map2color(Pseudotimes$Projected_Times_Cells, mypal)
    # paleta<-grDevices::colorRampPalette(colors=c("yellow", "orange", "red", "black"))
    # color<-paleta(200)
    plot(CellCoordinates[,1], CellCoordinates[,2], col=col1, pch=16, cex=1.5, ylab="Cooordinate 1", xlab="Coordinate 2")
  }
  else
  {
    print ("This function cannot plot pseudotime for datasets with more than 3 dimensions. If your dataset has more dimensions you may want to make a projection from those into 2 or 3 dimensions in order to visualize their pseudotimes")
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
#' @importFrom grDevices colorRampPalette
#' @importFrom fields image.plot
#' @importFrom stats as.dist cor hclust
plot_heatmaps_embedding <-function(Pseudotimes, EmbeddedTree, log_transform=F, cluster_genes=T)
{

  # # Testing
  # EmbeddedTree=EmbeddedTree
  # log_tranform=F
  # cluster_genes=T
  # # end testing

  CellByBranchTimes=c()
  CheckPseudotime=c()

  BranchesSequences=c()
  BranchesColors=c()

  selected_colors = c("forestgreen", "darkorchid", "darkorange3", "dodgerblue3", "firebrick3", "orange", "blue", "aquamarine", "wheat4", "brown", "gray")

  # Take cells from every branch and order them by pseudotime
  for(i in 1:length(Pseudotimes$Branches))
  {
    branch_i=Pseudotimes$Branches[[i]]
    aux=sort(Pseudotimes$Times_yk[branch_i], index.return=T)

    # generate a sequence of numbers to represent the limits between branches in the matrix
    if(i!=length(Pseudotimes$Branches))
    {
    BranchesSequences=c(BranchesSequences, rep(i, length(Pseudotimes$Branches[[i]])-1))
    BranchesColors=c(BranchesColors, rep(selected_colors[i], length(Pseudotimes$Branches[[i]])-1))
    }else{
      BranchesSequences=c(BranchesSequences, rep(i, length(Pseudotimes$Branches[[i]])))
      BranchesColors=c(BranchesColors, rep(selected_colors[i], length(Pseudotimes$Branches[[i]])))
    }
    CellByBranchTimes=c(CellByBranchTimes, branch_i[aux$ix])
  }

  CellByBranchTimes=unique(CellByBranchTimes)

  OrderedAverageNodes=EmbeddedTree$AveragedNodes[CellByBranchTimes,]
  OrderedInterpoaltedNodes=EmbeddedTree$Nodes[CellByBranchTimes,]

  if(log_transform==T)
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
    OrderedInterpoaltedNodes.hc <- stats::hclust(stats::as.dist(1-stats::cor(OrderedInterpoaltedNodes)))
    OrderedInterpoaltedNodes=OrderedInterpoaltedNodes[, OrderedInterpoaltedNodes.hc$labels[OrderedInterpoaltedNodes.hc$order]]
    OrderedAverageNodes=OrderedAverageNodes[, OrderedInterpoaltedNodes.hc$labels[OrderedInterpoaltedNodes.hc$order]]
  }


  my.colors = grDevices::colorRampPalette(c("black", "red","orange", "yellow", "white"))
  range_heatmap= c(min(c(range(OrderedAverageNodes)[1], range(OrderedInterpoaltedNodes)[1])), max(c(range(OrderedAverageNodes)[2], range(OrderedInterpoaltedNodes)[2])))
  fields::image.plot(t(OrderedAverageNodes), xlab="Genes", ylab="PseudoCells", col =my.colors(300) , zlim=range_heatmap, xaxt='n', yaxt='n', main="Averaged Gene Expression Profiles")
  fields::image.plot(t(OrderedInterpoaltedNodes), xlab="Genes", ylab="PseudoCells", col =my.colors(300) , zlim=range_heatmap, xaxt='n', yaxt='n', main="Imputed Gene Expression Profiles")

  EmbeddedMatrices= list(AveragedMatrix=OrderedAverageNodes, InterpolatedMatrix=OrderedInterpoaltedNodes)


  return(EmbeddedMatrices)
}


#' Plots Gene Expression on Cells Coordinates
#'
#' Plots the expression of a given gene on a given set of coordinates
#'
#' @param GeneName One of the elements in the Dataset$GeneNames object.
#' @param CellCoordinates Matrix with 2 or 3 coordinates per cell.
#' @param ExpressionMatrix Expression matrix
#' @param CoordLabels Which coordinate labels to use
#' @param Average wheter or not expression for cells will be averaged over the nearest knn neighbors. Default Average=FALSE
#' @param knn number of neighbors over which expression for the cells will be averaged in case Average=TRUE. Default: knn=5
#' @param Palette colors that will be used for the cells by default= Palette=colorRampPalette(colors=c("yellow", "orange", "red","black")))
#' @export
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats dist
#' @importFrom graphics plot
#' @importFrom rgl plot3d
plot_gene_on_map <- function(
  GeneName,
  CellCoordinates,
  ExpressionMatrix,
  CoordLabels=c("Component 1", "Component 2", "Component 3"),
  Average=F,
  knn=5,
  Palette=grDevices::colorRampPalette(colors=c("yellow", "orange", "red","black")))
{
  if (!dim(CellCoordinates)[2] %in% c(2,3)) {
    stop("Dimensions of CellCoordinates need to be 2 or 3.")
  }

  paleta = Palette
  color = paleta(400)
  Gene = which(colnames(ExpressionMatrix)==GeneName)
  ExpressionGene = as.numeric(ExpressionMatrix[,Gene])

  # In order to color the range properly, all values are applied an offset to be equal or greater to 0
  if(min(ExpressionGene)<0){ExpressionGene=ExpressionGene+min(ExpressionGene*-1)}

  if(Average==T)
  {
    ExpressionGeneAveraged=c()
    knn=5
    Distances=as.matrix(stats::dist(CellCoordinates, method = "euclidean", diag = FALSE, upper = TRUE, p = 2))

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

#' Plots a flattened 2 dimensional version of the Elastic Tree Nodes.
#'
#' Given an Elastic Tree it plots a 2 dimensional version of the Elastic Tree. If given annotation for each cell,
#' each node of the elastic tree will be a pie chart of the different cell identities that were assigned to the node.
#'
#' @param ElasticTree One of the elements in the Dataset$GeneNames object.
#' @param cell_annot Cell annotation to show in the nodes of the elastic tree. Default is the MERLoT branch identity.
#' @param layout Layout for the visualization of the elastic tree. Default is to use the Kamada-Kawai algorithm.
#' @param legend_position The position of the legend.
#' @param node_size If annotation was provided, set if the size of each pie chart will be uniform ("topology") or depend on the number of cells assigned to the node ("cells")
#' @param node_radius Sets the default node radius if \code{node_size = "topology"}
#' @param cols A vector of colors to use for the annotation.
#' @export
#'
#' @importFrom igraph graph_from_adjacency_matrix layout_with_kk
#' @importFrom graphics plot legend
#' @importFrom LSD distinctcolors
plot_flattened_tree <- function(ElasticTree,
                                cell_annot=NULL,
                                layout=NULL,
                                legend_position="topright",
                                node_size=c("topology", "cells"),
                                node_radius=12,
                                cols=NULL)
{
  node_size <- match.arg(node_size)
  NumberOfNodes <- dim(ElasticTree$Nodes)[1]
  N <- length(ElasticTree$Topology$Endpoints)
  B <- length(ElasticTree$Topology$Branchpoints)
  
  if (is.null(cols)) {
    selected_colors <- c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid",
                         "darkorange3", "orange", "blue", "aquamarine", "magenta",
                         "brown", "gray", "wheat1", "azure4", "lightsalmon4",
                         "navy", "sienna1", "gold4", "red4", "violetred", "black")
  } else {
    selected_colors <- cols
  }
  
  # build an adjacency matrix of the nodes and translate it to a graph
  graph_yk <- get_node_graph(ElasticTree)

  # calculate a viewable layout, if not given
  if (is.null(layout)) {
    l <- igraph::layout_with_kk(graph_yk, dim = 2)
  } else {
    l <- layout
  }

  # first decide what plot we are making
  # if cell_annot is empty, we are making the normal plot
  # if cell_annot has values we are making pie charts
  # if we are making pie charts, their size is either uniform
  # or it depends on how many cells belong to each pie
  if (is.null(cell_annot)) {
    nodes_colors=rep("black", dim(ElasticTree$Nodes)[1])

    for(i in 1:length(ElasticTree$Branches)) {
      nodes_colors[ElasticTree$Branches[[i]]]=selected_colors[i]
    }
    nodes_colors[1:N]="red"
    nodes_colors[(N+1):(N+B)]="skyblue"
    # nodes_colors= c(
    #   rep("red", length(ElasticTree$Topology$Endpoints)),
    #   rep("skyblue", length(ElasticTree$Topology$Branchpoints)),
    #   rep("black", NumberOfNodes -length(ElasticTree$Topology$Endpoints)-length(ElasticTree$Topology$Branchpoints))
    # )

    nodes_labels=c(seq(1, N + B), rep(NA, NumberOfNodes - N - B))
    nodes_sizes=c(rep(6, N + B), rep(4, NumberOfNodes - N - B))
    legend_labs <- c("Endpoints",
                     "Branchpoints",
                     paste("Branch ", 1:length(ElasticTree$Branches)))
    legend_cols <- c("red",
                     "skyblue",
                     selected_colors[1:length(ElasticTree$Branches)])
    
    plot(graph_yk, layout=l,
         vertex.label=nodes_labels,
         vertex.size=nodes_sizes,
         vertex.color=nodes_colors,
         vertex.label.cex=0.6)
    legend(x=legend_position,
           legend=legend_labs,
           col=legend_cols,
           pch=16)
  } else {
    annot_levels <- levels(as.factor(cell_annot))
    if (length(annot_levels) <= length(selected_colors)) {
      pie_cols <- c("white", selected_colors[1:length(annot_levels)])
    } else {
      print(paste("Supplied colors not enough; creating LSD palette with", length(annot_levels), "colors"))
      tmp <- LSD::distinctcolors(length(annot_levels), show = FALSE, bw = TRUE)
      pie_cols <- c("white", tmp)
    }
    annot_levels <- c("empty", annot_levels)
    pie_values <- lapply(1:NumberOfNodes, function(x) table(annot_levels) - 2)
    pie_cols <- pie_cols[order(annot_levels)]
    annot_levels <- annot_levels[order(annot_levels)]

    for (n in 1:NumberOfNodes) {
      cells_n <- (ElasticTree$Cells2TreeNodes[, 2] == n)
      for (real in annot_levels) {
        pie_values[[n]][real] <- sum(cell_annot[cells_n] == real)
      }
      if (sum(pie_values[[n]]) == 0) {
        pie_values[[n]]["empty"] <- 1
      }
    }

    nodes_labels=c(seq(1, N + B), rep(NA, NumberOfNodes - N - B))

    if (node_size == "topology") {
      nodes_sizes=c(rep(node_radius + 3, N + B), rep(node_radius, NumberOfNodes - N - B))
    } else if (node_size == "cells") {
      avg_cells_per_node <- length(cell_annot) / NumberOfNodes
      nodes_sizes = rep(5, NumberOfNodes)
      cells <- rep(0, NumberOfNodes)
      for (n in 1:NumberOfNodes) {
        cells[n] <- sum((ElasticTree$Cells2TreeNodes[, 2] == n))
        nodes_sizes[n] <- nodes_sizes[n] * (cells[n] / avg_cells_per_node)
      }
      sizes <- seq(0, max(cells), 10)[-1]
      size_length <- length(sizes)
      keep <- c(1,
                as.integer(size_length/4),
                as.integer(size_length/2),
                as.integer(size_length) * 0.75,
                size_length)
    }

    plot(graph_yk, layout=l,
         vertex.label=nodes_labels,
         vertex.size=nodes_sizes,
         vertex.shape="pie",
         vertex.pie=pie_values,
         vertex.pie.color=list(pie_cols),
         vertex.label.cex=0.6)
    l <- legend(x=legend_position,
                legend=annot_levels,
                col="black",
                pt.bg=pie_cols,
                pch=21)
    if (node_size == "cells") {
      a <- legend(x = l$rect$left, y=l$rect$top - l$rect$h,
                  legend=sizes[keep],
                  pt.cex= keep,
                  col='black')
      x <- (a$text$x + a$rect$left) / 2
      y <- a$text$y
      graphics::symbols(x,y,circles=sizes[keep] * 5 / (200 * avg_cells_per_node) ,inches=FALSE,add=TRUE,bg='orange')
    }
  }
  # legend(x="bottomright", inset=c(-0.16,0), legend=paste("Branch ", 1:length(ElasticTree$Branches)), col = selected_colors[1:length(ElasticTree$Branches)], pch=c(16))
}

#' Calculates an adjacency matrix for the support nodes of an elastic tree and returns the corresponding graph.
#'
#' @param ElasticTree the tree to convert
#'
#' @return igraph graph
#'
#' @export
get_node_graph <- function(ElasticTree) {
  NumberOfNodes <- dim(ElasticTree$Nodes)[1]
  EdgesTree <- matrix(0, NumberOfNodes, NumberOfNodes)
  for(i in (1: dim(ElasticTree$Edges)[1])) {
    EdgesTree[ElasticTree$Edges[i, 1], ElasticTree$Edges[i, 2]] <- 1
  }
  graph_yk <- igraph::graph_from_adjacency_matrix(EdgesTree, mode = "undirected")
  return(graph_yk)
}


#' Plots the gene expression level of a gene on a flattened 2 dimensional version of the Elastic Tree Nodes
#'
#' Given an Elastic Tree it plots a 2 dimensional version of the Elastic Tree.
#'
#' @param ElasticTree One of the elements in the Dataset$GeneNames object.
#' @param GeneName Name of the gene to be visualized. The gene name is one of the colnames in the EmbeddedTree$Nodes object
#' @param legend_position The position of the legend.
#'
#' @export
#'

plot_flattened_tree_gene_expression <- function(ElasticTree, GeneName, legend_position="topright")
{
  NumberOfNodes=dim(ElasticTree$Nodes)[1]

  selected_colors = c("forestgreen", "firebrick3", "dodgerblue3", "darkorchid", "darkorange3", "orange", "blue", "aquamarine", "magenta", "brown", "gray", "wheat1", "azure4", "lightsalmon4", "navy", "sienna1", "gold4", "red4", "violetred")

  EdgesTree=matrix(0, dim(ElasticTree$Nodes)[1], dim(ElasticTree$Nodes)[1])

  for(i in (1: dim(ElasticTree$Edges)[1]))
  {
    EdgesTree[ElasticTree$Edges[i, 1], ElasticTree$Edges[i, 2]]=1
  }

  graph_yk=graph_from_adjacency_matrix(EdgesTree, mode = "undirected")
  l=layout_with_kk(graph_yk, dim = 2)

  nodes_colors=rep("black", dim(ElasticTree$Nodes)[1])

  for(i in 1:length(ElasticTree$Branches))
  {
    nodes_colors[ElasticTree$Branches[[i]]]=selected_colors[i]
  }

  nodes_colors[1:length(ElasticTree$Topology$Endpoints)]="red"
  nodes_colors[(length(ElasticTree$Topology$Endpoints)+1):(length(ElasticTree$Topology$Endpoints)+length(ElasticTree$Topology$Branchpoints))]="skyblue"

  nodes_labels=c(seq(1, length(ElasticTree$Topology$Endpoints) + length(ElasticTree$Topology$Branchpoints)), rep(NA, NumberOfNodes -length(ElasticTree$Topology$Endpoints)-length(ElasticTree$Topology$Branchpoints)))
  nodes_sizes=c(rep(6, length(ElasticTree$Topology$Endpoints) + length(ElasticTree$Topology$Branchpoints)), rep(4, NumberOfNodes -length(ElasticTree$Topology$Endpoints)-length(ElasticTree$Topology$Branchpoints)))

  colnames(ElasticTree$Nodes)
  Gene = ElasticTree$Nodes[, GeneName]
  mypal <- colorRampPalette( c( "yellow", "orange", "darkorange", "red", "black" ) )( length(Gene) )
  nodes_colors=map2color(Gene, mypal)
  plot(graph_yk, layout=l, vertex.label=nodes_labels, vertex.size=nodes_sizes, vertex.color=nodes_colors, vertex.label.cex=0.6, main=paste("gene expression: ", GeneName))
}

