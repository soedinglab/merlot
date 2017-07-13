#' Calculate Elastic Tree
#'
#'Calculates an elastic tree using a scaffold tree to initialize it. A set of N_yk nodes are included one byone into a tree structure that minimizes an error and an energetic functions to the cell coordinates. The initialization with the scaffold tree ensures that the correct number of endpoints and branchpoints and their connectivity are preserved in the elastic tree interpolation.
#' @param Scaffoldtree scaffoldTree calculated by the CalculateScaffoldTree function.
#' @param N_yk number of nodes for the elastic principal tree
#' @param lambda_0 principal elastic tree energy function parameter.
#' @param mu_0 principal elastic tree energy function parameter.
#' @return ElasticTree
#' @export
CalculateElasticTree <- function(ScaffoldTree, N_yk=100, input="topology", lambda_0=2.03e-09, mu_0=0.00625, FixEndpoints=F, plot=F)
{
  # Testing
  # Default parameters taken from adjustment in real datasets with 100 N_yks
  # lambda_0=2.03e-09
  # mu_0=0.00625
  #  End testing

  Coords=c()
  Edges=c()

  # scaling according to N_yk
  mu=(N_yk-1)*mu_0
  lambda=((N_yk-2)**3)*lambda_0

  # I apply unique in case there are repeated nodes, which is a consequence of having trifurcations or higher order connections
  TopologyNodes=unique(c(ScaffoldTree$Endpoints, ScaffoldTree$Branchpoints))
  # Calculate connectivity in between branches
  TopologyEdges=c()
  TopologyCoords=ScaffoldTree$CellCoordinates[TopologyNodes,]

  for(i in (1:dim(ScaffoldTree$Branches)[1]))
  {
    TopologyEdges=rbind(TopologyEdges, c(which(TopologyNodes==ScaffoldTree$Branches[i,1])[1], which(TopologyNodes==ScaffoldTree$Branches[i,2])[1]))
  }

  Coords=TopologyCoords
  Edges=TopologyEdges

  # Given initial topology, calculate connectivity in between branches
  TopologyNodes=c(ScaffoldTree$Endpoints, ScaffoldTree$Branchpoints)
  EmbeddedTopology=c()
  for(i in (1:dim(ScaffoldTree$Branches)[1]))
  {
    EmbeddedTopology=rbind(EmbeddedTopology, (c(which(TopologyNodes==ScaffoldTree$Branches[i,1])[1], which(TopologyNodes==ScaffoldTree$Branches[i,2])[1])))
  }
  Edges=EmbeddedTopology
  ElasticTree <- computeElasticPrincipalGraph(Data = ScaffoldTree$CellCoordinates, NumNodes = N_yk,
                                              NodesPositions = Coords, Edges = Edges,
                                              Method = 'CurveConfiguration', EP=lambda, RP=mu)

  # Unlist the ElasticTree structure
  ElasticTree=ElasticTree[[1]]

  # Add the topology in the tree structure
  ETreeTopology= list(Endpoints=seq(from=1, to=length(ScaffoldTree$Endpoints), by=1), Branchpoints=seq(from=(length(ScaffoldTree$Endpoints) +1), to= (length(ScaffoldTree$Endpoints)) + length(ScaffoldTree$Branchpoints), by=1))
  ElasticTree=c(ElasticTree, Topology=1)
  ElasticTree$Topology <- ETreeTopology

  # Add Input Coordinates to the tree structure
  ElasticTree=c(ElasticTree, CellCoords=1)
  ElasticTree$CellCoords <- ScaffoldTree$CellCoordinates

  # Add Tree Connectivity
  ElasticTree=c(ElasticTree, Connectivity=1)
  ElasticTree$Connectivity  <- Edges

  # Assign nodes from the embedded tree to the different branches and save the structure in the tree structure
  BranchesNodes=list()
  for (i in 1: dim(EmbeddedTopology)[1])
  {
    path=EmbeddedTopology[i,1]
    first=EmbeddedTopology[i,1]
    while (first!=EmbeddedTopology[i,2])
    {
      edge=ElasticTree$Edges[which(ElasticTree$Edges[,1]==first),]
      second=edge[2]
      path=c(path, second)
      first=second
    }
    BranchesNodes[[i]] <- path
  }

  cell2yk=c()
  for (i in 1:dim(ScaffoldTree$CellCoordinates)[1])
  {
    cell_i=matrix(ScaffoldTree$CellCoordinates[i,], nrow=1)
    dist_cell_i=as.matrix(dist(rbind(cell_i, ElasticTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk=rbind(cell2yk, c(i, closest_yk))
  }

  # Add cells 2 yks mapping
  ElasticTree=c(ElasticTree, Cells2TreeNodes=1)
  ElasticTree$Cells2TreeNodes <- cell2yk

  # assign cells to the branches
  allnodes=c()
  for(i in 1:length(BranchesNodes))
  {
    allnodes=c(allnodes, BranchesNodes[[i]])
  }
  allnodes=sort(unique(allnodes))
  cells_branchs_assigments=c()
  for (i in 1:length(BranchesNodes))
  {
    cells_branchs_assigments[which(cell2yk[,2] %in% BranchesNodes[[i]])]=i
  }

  ElasticTree=c(ElasticTree, Branches=1)
  ElasticTree$Branches <- BranchesNodes

  ElasticTree=c(ElasticTree, Cells2Branches=1)
  ElasticTree$Cells2Branches <- cells_branchs_assigments

  # Fixing endpoints coordinates
  if(FixEndpoints==T)
  {
    ElasticTree$Nodes[1:length(ScaffoldTree$Endpoints),]=ScaffoldTree$CellCoordinates[ScaffoldTree$Endpoints,]
  }

  return (ElasticTree)
}

#' Gene Space Embedding
#'
#' Embeds an elastic principal tree on the gene expression space from which the initial cell coordinates were calculated. The function receives the original Expression Matrix from which a lowe dimensional manifold was calculated and the coordinates from the low-dimensional manifold.
#'
#' @param ExpressionMatrix Gene expression matrix on which the tree will be embedded
#' @param ElasticTree elastic tree to be embedded on the gene expression space
#' @param increaseFactor factor by which the principal elastic tree energy function parameters will be increased for the embedding. By default this number is 10.
#' @export

GenesSpaceEmbedding <- function(ExpressionMatrix, ElasticTree,  lambda_0=2.03e-09, mu_0=0.00625, increaseFactor=10)
{
  # The number of nodes for the embedding tree is the same as the ones for the input low dimensional one
  N_yk=dim(ElasticTree$Nodes)[1]

  # scaling according to N_yk and multiply the values by the increase factor
  mu=(N_yk-1) * mu_0 * increaseFactor
  lambda=((N_yk-2)**3) * lambda_0 * increaseFactor

  # Map each xn (cells) to the closest y_k (tree node)
  cell2yk=c()
  for (i in 1:dim(ExpressionMatrix)[1])
  {
    cell_i=matrix(ElasticTree$CellCoords[i,], nrow = 1)
    dist_cell_i=as.matrix(dist(rbind(cell_i, ElasticTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk=rbind(cell2yk, c(i, closest_yk))
  }

  # count number of xn per y_k
  yk_counts=c()
  yk_profiles=c()

  # calculate the transcriptional profile for y_ks
  for (i in 1:dim(ElasticTree$Nodes)[1])
  {
    # count cells associated to y_ki
    yk_counts=c(yk_counts, length(cell2yk[which(cell2yk[,2]==i),]))

    # when more than 1 x_n is mapped to y_ki
    if(is.matrix(ExpressionMatrix[which(cell2yk[,2]==i),]))
    {
      yk_profiles=rbind(yk_profiles, colMeans(ExpressionMatrix[which(cell2yk[,2]==i),]))
    }
    # if only 1 or 0 x_n is mapped to y_ki
    else  if(length(ExpressionMatrix[which(cell2yk[,2]==i),])==dim(ExpressionMatrix)[2])
    {
      yk_profiles=rbind(yk_profiles, ExpressionMatrix[which(cell2yk[,2]==i),])
    }
  }

  # set NA values for y_k with no mapped x_n to 0
  yk_profiles=as.matrix(yk_profiles)
  yk_profiles[which(is.na(yk_profiles))]=0
  length(which(yk_counts==0))

  # # # Work in progress
  # for(i in 1:N_yk)
  # {
  #   if(length(yk_profiles[i,is.na(yk_profiles[i,])])==length(yk_profiles[i,]))
  #   {
  #     # search previous non null node
  #     previous_node=i
  #     while (length(yk_profiles[previous_node,is.na(yk_profiles[previous_node,])])==length(yk_profiles[i,]))
  #     {
  #       previous_node=previous_node-1
  #     }
  #
  #     # search previous non null node
  #     next_node=i
  #     while (length(yk_profiles[next_node,is.na(yk_profiles[next_node,])])==length(yk_profiles[i,]))
  #     {
  #       next_node=next_node+1
  #     }
  #
  #     yk_profiles[i,]=(yk_profiles[previous_node,] + yk_profiles[next_node,])/2
  #   }
  # }
  # # # end work in progress

  # calculate elastic tree for Expression Data using Y_ks profiles and edges
  EmbeddedTree <- computeElasticPrincipalGraph(Data = data.matrix(ExpressionMatrix), NumNodes = N_yk,
                                               NodesPositions = yk_profiles, Edges = ElasticTree$Edges,
                                               Method = 'CurveConfiguration',  EP=lambda, RP=mu)

  # Unlist EmbeddedTree structure
  EmbeddedTree=EmbeddedTree[[1]]

  # reassigning cells to nodes on the full dimensional space
  cell2yk_post=c()
  for (i in 1:dim(ExpressionMatrix)[1])
  {
    cell_i= matrix(ExpressionMatrix[i,], nrow=1)
    dist_cell_i=as.matrix(dist(rbind(cell_i, EmbeddedTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk_post=rbind(cell2yk_post, c(i, closest_yk))
  }

  allnodes=1:N_yk
  cells_branchs_assigments=c()
  for (i in 1:length(ElasticTree$Branches))
  {
    cells_branchs_assigments[which(cell2yk_post[,2] %in% ElasticTree$Branches[[i]])]=i
  }

  # Add tree topology to the tree structure
  EmbeddedTree=c(EmbeddedTree, Topology=1)
  EmbeddedTree$Topology <- ElasticTree$Topology

  # Add tree connectivity to the tree structure
  EmbeddedTree=c(EmbeddedTree, Connectivity=1)
  EmbeddedTree$Connectivity <- ElasticTree$Connectivity

  # Add gene names to columns
  colnames(EmbeddedTree$Nodes) <-colnames(ExpressionMatrix)

  # Add Input Coordinates to the tree structure
  EmbeddedTree=c(EmbeddedTree, CellCoords=1)
  EmbeddedTree$CellCoords <- ExpressionMatrix

  # Add cells 2 yks mapping
  EmbeddedTree=c(EmbeddedTree, Cells2TreeNodes=1)
  EmbeddedTree$Cells2TreeNodes <- cell2yk_post

  # Add cells 2 branches mapping
  EmbeddedTree=c(EmbeddedTree, Cells2Branches=1)
  EmbeddedTree$Cells2Branches <- cells_branchs_assigments

  # Assigning branches from input tree to this one
  EmbeddedTree=c(EmbeddedTree, Branches=1)
  EmbeddedTree$Branches <- ElasticTree$Branches

  # Calculate diffusion map for resulting tree taking out vectors equal to 0
  yk_profiles_nozeros=yk_profiles[which(!(seq(1,N_yk,1) %in% which(rowMeans((yk_profiles))==0))),]

  # testing
  EmbeddedTree=c(EmbeddedTree, AveragedNodes=1)
  EmbeddedTree$AveragedNodes <- yk_profiles

  # dif_yk_profiles <- DiffusionMap(yk_profiles_nozeros, verbose = T)
  dif_elastic_yk <- DiffusionMap(EmbeddedTree$Nodes,  density.norm = T)

  return (EmbeddedTree)
}
