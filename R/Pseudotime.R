#' Calculate Pseudotime
#'
#' Calculates pseudotime for the cells either on the low dimensional space or the full gene expression space. Users can define which is the endpoint or cell that will be used as the initial time (t0) respect to which the pseudotimes to all cells will be calculated.
#'
#' @param InputTree elastic tree according to which the pseudotime will be calculated. It can be a low dimensional elastic tree as the one calculated by CalculateElasticTree() or a high dimensional elastic tree like the one calculated by GeneSpaceEmbedding()
#' @param T0 which endpoint will be used at time zero
#' @param C0 which cell will be used as time zero
#' @export
#'
#' @importFrom stats dist
#' @importFrom igraph graph_from_adjacency_matrix shortest.paths
CalculatePseudotimes <- function (InputTree, plot=F, plotdim, T0=1, C0=NULL)
{
  # Testing
  # InputTree=ElasticTree3DTopology
  # End Testing

  N_yk=dim(InputTree$Nodes)[1]
  # Calculate adjacency matrix for the embedded tree
  Distances_yk=matrix(0, N_yk, N_yk)

  for(i in (1: dim(InputTree$Edges)[1]))
  {
    Distances_yk[InputTree$Edges[i, 1], InputTree$Edges[i, 2]]=1
    Distances_yk[InputTree$Edges[i, 2], InputTree$Edges[i, 1]]=1
  }

  Graph_yk=igraph::graph_from_adjacency_matrix(Distances_yk)
  Times_yk=igraph::shortest.paths(Graph_yk)
  # shortest path distances from the 1st Endpoint to all the other cells. The first endpoint is always the first element in the array
  if(!is.null(C0))
  {
    T0=InputTree$Cells2TreeNodes[which(InputTree$Cells2TreeNodes[,1]==C0), 2]
  }

  Times_yk=Times_yk[T0,]

  # Map cells to the closest yk node and assign its pseudotime
  cell2yk=c()
  Times_cells=c()
  Proyected_Times_Cells=c()

  for (i in 1:dim(InputTree$CellCoords)[1])
  {
    cell_i=matrix(InputTree$CellCoords[i,], nrow=1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, InputTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))

    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    second_closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[3]-1

    # calculate extra pseudotime in the grid
    # x= (c^2+a^2+b^2 )/ 2*a
    b=sort(dist_cell_i[,1], index.return=T)$x[2]
    c=sort(dist_cell_i[,1], index.return=T)$x[3]
    a=stats::dist(rbind(matrix(InputTree$Nodes[closest_yk,], nrow=1), matrix(InputTree$Nodes[second_closest_yk,], nrow=1)),  method = "euclidean", diag=F, p=2)
    x= (a^2+b^2-c^2 )/ (2*a)
    extra_time=x/a

    cell2yk=rbind(cell2yk, c(i, closest_yk))
    Times_cells=c(Times_cells, Times_yk[closest_yk])
    Proyected_Times_Cells=c(Proyected_Times_Cells, (Times_yk[closest_yk] + extra_time))
  }
  # delete names from vector
  names(Proyected_Times_Cells) = NULL

  allnodes=1:N_yk
  cells_branchs_assigments=c()
  for (i in 1:length(InputTree$Branches))
  {
    cells_branchs_assigments[which(InputTree$Cells2TreeNodes[,2] %in% InputTree$Branches[[i]])]=i
  }

  #Assign expression from closest y_kn node only valid if an Embedded Tree is used as input
  if(dim(InputTree$CellCoords)[2]==dim(InputTree$Nodes)[2])
  {
    PseudoExpressionMatrix=InputTree$CellCoords
    for(i in 1:dim(InputTree$CellCoords)[1])
    {
      cell2yk[which(cell2yk[,1]==i),2]
      PseudoExpressionMatrix[i,]=InputTree$Nodes[cell2yk[which(cell2yk[,1]==i),2],]
    }
    Pseudotimes= list(Times_yk=Times_yk, Times_cells=Times_cells, Proyected_Times_Cells=Proyected_Times_Cells, Branches=InputTree$Branches, Cells2TreeNodes=cell2yk, Cells2Branches=cells_branchs_assigments, PseudoExpressionMatrix=PseudoExpressionMatrix, T0=T0, C0=C0)
  }
  else
  {
    Pseudotimes= list(Times_yk=Times_yk, Times_cells=Times_cells, Proyected_Times_Cells=Proyected_Times_Cells, Branches=InputTree$Branches, Cells2TreeNodes=cell2yk, Cells2Branches=cells_branchs_assigments, T0=T0, C0=C0)
  }
  return (Pseudotimes)
}
