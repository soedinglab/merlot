#' Calculate Elastic Tree
#'
#'Calculates an elastic tree using a scaffold tree to initialize it. A set of N_yk nodes are included one byone into a tree structure that minimizes an error and an energetic functions to the cell coordinates. The initialization with the scaffold tree ensures that the correct number of endpoints and branchpoints and their connectivity are preserved in the elastic tree interpolation.
#' @param ScaffoldTree scaffoldTree calculated by the CalculateScaffoldTree function.
#' @param N_yk number of nodes for the elastic principal tree
#' @param lambda_0 principal elastic tree energy function parameter.
#' @param mu_0 principal elastic tree energy function parameter.
#' @param FixEndpoints Whether or not the end points coordinates are fixed
#' @param plot Whether or not plots are produced
#' @param NBranchScaffoldNodes Whether or not to add middle scaffold nodes
#' @param NCores number of cpu cores to be used for the calculation
#' @return ElasticTree
#' @export
#'
#' @importFrom stats dist
#' @importFrom ElPiGraph.R computeElasticPrincipalCurve
#' @importFrom igraph graph_from_adjacency_matrix get.shortest.paths
CalculateElasticTree <- function(ScaffoldTree, N_yk=100, lambda_0=0.80e-09, mu_0=0.00250, FixEndpoints=F, plot=F, NBranchScaffoldNodes = 1, NCores=1)
{
  # Testing
  # Default parameters taken from adjustment in real datasets with 100 N_yks
  # N_yk=100
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
  TopologyCoords=ScaffoldTree$CellCoordinates[TopologyNodes,]
  # Given initial topology, calculate connectivity in between branches
  TopologyNodes=c(ScaffoldTree$Endpoints, ScaffoldTree$Branchpoints)
  TopologyEdges=c()

  if(length(ScaffoldTree$Endpoints)==2)
  {
    TopologyEdges=rbind(TopologyEdges, (c(which(TopologyNodes==ScaffoldTree$Branches[1])[1], which(TopologyNodes==ScaffoldTree$Branches[2])[1])))
  }else
  {
    for(i in (1:dim(ScaffoldTree$Branches)[1]))
    {
      TopologyEdges=rbind(TopologyEdges, (c(which(TopologyNodes==ScaffoldTree$Branches[i,1])[1], which(TopologyNodes==ScaffoldTree$Branches[i,2])[1])))
    }
  }


  TopologyCoordsAux=TopologyCoords
  TopologyEdgesAux=TopologyEdges
  TopologyNodesAux=TopologyNodes

  # ---------------------------------------
  # -------Add branch middle scaffold nodes
  # ---------------------------------------
  if(NBranchScaffoldNodes)
  {
    TopologyEdgesAux=c()

    for(i in 1:dim(ScaffoldTree$Branches)[1])
    {
      # Distance between to branch extreme cells
      extremepointsDistance=ScaffoldTree$DijkstraDistances[ScaffoldTree$Branches[i,1], ScaffoldTree$Branches[i,2]]
      # Cells between those two branches
      nodesBranchi=ScaffoldTree$SkeletonNodes[which(ScaffoldTree$SkeletonNodes %in% ScaffoldTree$Cells2BranchesAssignments[[i]])]
      # take the endpoints out of the array
      nodesBranchi=nodesBranchi[which(!(nodesBranchi %in% ScaffoldTree$Branches[i,]))]

      minDistanceCenter=min((ScaffoldTree$DijkstraDistances[ScaffoldTree$Branches[i,1],nodesBranchi])-(extremepointsDistance/2))
      Distances2Nodes=ScaffoldTree$DijkstraDistances[ScaffoldTree$Branches[i,1],nodesBranchi]

      # New ScaffoldNode
      NewScaffoldNode=which(abs(Distances2Nodes-(extremepointsDistance/2)) == min(abs(Distances2Nodes-extremepointsDistance/2)))
      TopologyCoordsAux=rbind(TopologyCoordsAux, ScaffoldTree$CellCoordinates[NewScaffoldNode,])
      NewScaffoldNodePosition=length(TopologyNodesAux)+1
      TopologyNodesAux=c(TopologyNodesAux, NewScaffoldNode)
      # NedEdges
      Edge1=c(which(TopologyNodesAux==ScaffoldTree$Branches[i,1]), NewScaffoldNodePosition)
      Edge2=c(which(TopologyNodesAux==ScaffoldTree$Branches[i,2]), NewScaffoldNodePosition)
      TopologyEdgesAux=rbind(TopologyEdgesAux, Edge1)
      TopologyEdgesAux=rbind(TopologyEdgesAux, Edge2)
      rownames(TopologyEdgesAux)<-NULL
    }
  }

  ElasticTree=ElPiGraph.R::computeElasticPrincipalCurve(
    X = ScaffoldTree$CellCoordinates, NumNodes = N_yk,
    InitNodePositions = TopologyCoordsAux, InitEdges = TopologyEdgesAux,
    Lambda = lambda, Mu = mu, Do_PCA = F, verbose = F, drawAccuracyComplexity = F,
    drawPCAView = F, drawEnergy = F, Mode = 1, n.cores = NCores)

  # Unlist the ElasticTree structure
  ElasticTree=ElasticTree[[1]]
  names(ElasticTree)[1]="Nodes"
  ElasticTree$Edges=ElasticTree$Edges$Edges

  # Add the topology in the tree structure
  ETreeTopology= list(Endpoints=seq(from=1, to=length(ScaffoldTree$Endpoints), by=1), Branchpoints=seq(from=(length(ScaffoldTree$Endpoints) +1), to= (length(ScaffoldTree$Endpoints)) + length(ScaffoldTree$Branchpoints), by=1))
  ElasticTree=c(ElasticTree, Topology=1)
  ElasticTree$Topology <- ETreeTopology

  # Coords for the topology elements
  ElasticTree=c(ElasticTree, EndpointsCoords=1)
  ElasticTree$EndpointsCoords=ElasticTree$Nodes[ElasticTree$Topology$Endpoints,]
  ElasticTree=c(ElasticTree, BranchpointsCoords=1)
  ElasticTree$BranchpointsCoords=ElasticTree$Nodes[ElasticTree$Topology$Branchpoints,]

  # Add Input Coordinates to the tree structure
  ElasticTree=c(ElasticTree, CellCoords=1)
  ElasticTree$CellCoords <- ScaffoldTree$CellCoordinates

  # Add Tree Connectivity
  ElasticTree=c(ElasticTree, Connectivity=1)
  ElasticTree$Connectivity  <- Edges

  # Assign nodes from the embedded tree to the different branches and save the structure in the tree structure
  BranchesNodes=list()
  EdgesTree=matrix(0, N_yk, N_yk)

  for(i in (1: dim(ElasticTree$Edges)[1]))
  {
    EdgesTree[ElasticTree$Edges[i, 1], ElasticTree$Edges[i, 2]]=1
    EdgesTree[ElasticTree$Edges[i, 2], ElasticTree$Edges[i, 1]]=1
  }
  Graph_yk=igraph::graph_from_adjacency_matrix(EdgesTree)
  for(i in 1:dim(TopologyEdges)[1])
  {
    path_brach_i=igraph::get.shortest.paths(Graph_yk,from = TopologyEdges[i,1], to = TopologyEdges[i,2])
    BranchesNodes[[i]]=path_brach_i$vpath[[1]]
  }

  cell2yk=c()

  if(length(ScaffoldTree$Endpoints)==2)
  {
    cell_i=matrix(ScaffoldTree$CellCoordinates[1,], nrow=1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, ElasticTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk=rbind(cell2yk, c(1, closest_yk))
  }else
  {
    for (i in 1:dim(ScaffoldTree$CellCoordinates)[1])
    {
      cell_i=matrix(ScaffoldTree$CellCoordinates[i,], nrow=1)
      dist_cell_i=as.matrix(stats::dist(rbind(cell_i, ElasticTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
      #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
      closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
      cell2yk=rbind(cell2yk, c(i, closest_yk))
    }
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

#' Duplicate Elastic Tree Nodes
#' Introduces intermediate nodes in between nodes being part of an edge in the ElasticTree structure
#'
#' @param ElasticTree Elastic Tree to which nodes will be added
#'
#' @return ElasticTree
#' @export
#'
#' @importFrom stats dist
DuplicateTreeNodes <- function(ElasticTree)
{
  ElasticTree2=ElasticTree
  N_yk=dim(ElasticTree$Nodes)[1]
  NewEdges=c()
  NewNodes=c()

  for(i in 1:dim(ElasticTree$Edges)[1])
    # for(i in 1:3)
  {
    Edge_i=ElasticTree$Edges[i,]
    coords_node1=ElasticTree$Nodes[Edge_i[1],]
    coords_node2=ElasticTree$Nodes[Edge_i[2],]
    coords_interpolate=coords_node1 + ((coords_node2 - coords_node1)/2)

    # node N_yk + i
    NewEdges=rbind(NewEdges, c(Edge_i[1], N_yk+i))
    NewEdges=rbind(NewEdges, c(Edge_i[2], N_yk+i))
    NewNodes=rbind(NewNodes, coords_interpolate)


    # Find branch
    branch_edges=c()
    for(j in 1:length(ElasticTree2$Branches))
    {
      if(Edge_i[1] %in% ElasticTree2$Branches[[j]] && Edge_i[2] %in% ElasticTree2$Branches[[j]])
      {
        branch_edges=j
      }
    }

    if(which(ElasticTree2$Branches[[branch_edges]]==Edge_i[2]) < which(ElasticTree2$Branches[[branch_edges]]==Edge_i[1]))
    {
       Edge_i=rev(Edge_i)
    }

    ElasticTree2$Branches[[branch_edges]]=c(ElasticTree2$Branches[[branch_edges]][1:which(ElasticTree2$Branches[[branch_edges]]==Edge_i[1])], N_yk+i, ElasticTree2$Branches[[branch_edges]][which(ElasticTree2$Branches[[branch_edges]]==Edge_i[2]):length(ElasticTree2$Branches[[branch_edges]])])
  }

  ElasticTree2$Edges=NewEdges
  rownames(NewNodes)=c()
  ElasticTree2$Nodes=rbind(ElasticTree2$Nodes, NewNodes)

  # Reassing cells to nodes
  # reassigning cells to nodes on the full dimensional space
  cell2yk_post=c()
  for (i in 1:dim(ElasticTree2$CellCoords)[1])
  {
    cell_i= matrix(ElasticTree2$CellCoords[i,], nrow=1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, ElasticTree2$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk_post=rbind(cell2yk_post, c(i, closest_yk))
  }

  ElasticTree2$Cells2TreeNodes=cell2yk_post

  return(ElasticTree2)
}

#' Gene Space Embedding
#'
#' Embeds an elastic principal tree on the gene expression space from which the initial cell coordinates were calculated. The function receives the original Expression Matrix from which a lowe dimensional manifold was calculated and the coordinates from the low-dimensional manifold.
#'
#' @param ExpressionMatrix Gene expression matrix on which the tree will be embedded
#' @param ElasticTree elastic tree to be embedded on the gene expression space
#' @inheritParams CalculateElasticTree
#' @param increaseFactor_mu factor by which the principal elastic tree energy function parameters will be increased for the embedding.
#' @param increaseFactor_lambda factor by which the principal elastic tree energy function parameters will be increased for the embedding.
#'
#' @export
#'
#' @importFrom stats dist
#' @importFrom ElPiGraph.R computeElasticPrincipalCurve
GenesSpaceEmbedding <- function(ExpressionMatrix, ElasticTree,  lambda_0=2.03e-09, mu_0=0.00625, increaseFactor_mu=20, increaseFactor_lambda=20, NCores=1)
{
  # The number of nodes for the embedding tree is the same as the ones for the input low dimensional one
  N_yk=dim(ElasticTree$Nodes)[1]

  # scaling according to N_yk and multiply the values by the increase factor
  mu=(N_yk-1) * mu_0 * increaseFactor_mu
  lambda=((N_yk-2)**3) * lambda_0 * increaseFactor_lambda

  # Map each xn (cells) to the closest y_k (tree node)
  cell2yk=c()
  for (i in 1:dim(ExpressionMatrix)[1])
  {
    cell_i=matrix(ElasticTree$CellCoords[i,], nrow = 1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, ElasticTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk=rbind(cell2yk, c(i, closest_yk))
  }

  # count number of xn per y_k
  yk_counts <- rep(0, N_yk)
  yk_profiles <- matrix(data = 0, ncol = dim(ExpressionMatrix)[2], nrow = N_yk)
  colnames(yk_profiles) <- colnames(ExpressionMatrix)
  # calculate the transcriptional profile for y_ks
  for (i in 1:dim(ElasticTree$Nodes)[1])
  {
    # count cells associated to y_ki
    yk_counts[i] <- length(cell2yk[which(cell2yk[,2]==i),1])
    if (yk_counts[i] > 1) {
      yk_profiles[i,] <- colMeans(ExpressionMatrix[which(cell2yk[,2]==i),])
    } else if (yk_counts[i] == 1) {
      yk_profiles[i,] <- unlist(ExpressionMatrix[which(cell2yk[,2]==i),])
    }
  }

  # set NA values for y_k with no mapped x_n to 0
  yk_profiles[yk_counts == 0, ] <- 0

  CellCoordinates=data.matrix(ExpressionMatrix)
  InitialNodesCoordinates=yk_profiles
  InitialEdges= ElasticTree$Edges

  EmbeddedTree=ElPiGraph.R::computeElasticPrincipalCurve(
    X = CellCoordinates, NumNodes = N_yk,
    InitNodePositions = InitialNodesCoordinates, InitEdges = InitialEdges,
    Do_PCA = F, verbose = T, drawAccuracyComplexity = F,
    drawPCAView = F, drawEnergy = F, Lambda = lambda, Mu = mu, Mode = 1, n.cores = NCores)

  # Unlist EmbeddedTree structure
  EmbeddedTree=EmbeddedTree[[1]]
  names(EmbeddedTree)[1]="Nodes"
  EmbeddedTree$Edges=EmbeddedTree$Edges$Edges

  # reassigning cells to nodes in the full dimensional space
  cell2yk_post=c()
  for (i in 1:dim(ExpressionMatrix)[1])
  {
    cell_i= matrix(ExpressionMatrix[i,], nrow=1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, EmbeddedTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
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

  return (EmbeddedTree)
}


#' Calculate elastic tree with iterative constraints
#'
#'Calculates an elastic tree using a scaffold tree to initialize it. A set of N_yk nodes are included one by one into a tree structure that minimizes an error and an energetic functions to the cell coordinates. The initialization with the scaffold tree ensures that the correct number of endpoints and branchpoints and their connectivity are preserved in the elastic tree interpolation.
#' @param start_N_yk initial number of nodes with which the elastic tree will be calculated on top of which additional constraints will be added
#' @param step_N_yk number of nodes that will be interatively added to the elastic tree until the N_yk nodes are added. After each iteration the nodes will be added as additional constraints to the next iteration.
#' @inheritParams CalculateElasticTree
#' @return ElasticTree
#' @export
#'
#' @importFrom stats dist
#' @importFrom ElPiGraph.R computeElasticPrincipalGraph
CalculateElasticTreeConstrained <- function(ScaffoldTree, N_yk=150, start_N_yk=100, step_N_yk=50, lambda_0=2.03e-09, mu_0=0.00625, FixEndpoints=F, plot=F)
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
  TopologyCoords=ScaffoldTree$CellCoordinates[TopologyNodes,]
  # Given initial topology, calculate connectivity in between branches
  TopologyNodes=c(ScaffoldTree$Endpoints, ScaffoldTree$Branchpoints)
  TopologyEdges=c()

  if(length(ScaffoldTree$Endpoints)==2)
  {
    TopologyEdges=rbind(TopologyEdges, (c(which(TopologyNodes==ScaffoldTree$Branches[1])[1], which(TopologyNodes==ScaffoldTree$Branches[2])[1])))
  }else
  {
    for(i in (1:dim(ScaffoldTree$Branches)[1]))
    {
      TopologyEdges=rbind(TopologyEdges, (c(which(TopologyNodes==ScaffoldTree$Branches[i,1])[1], which(TopologyNodes==ScaffoldTree$Branches[i,2])[1])))
    }
  }

  N_yk_limits=(seq(from=start_N_yk, to=N_yk, by=step_N_yk))
  N_yk_limits=unique(c(N_yk_limits, N_yk))

  mu=(N_yk_limits[1]-1)*mu_0
  lambda=((N_yk_limits[1]-2)**3)*lambda_0

  ElasticTree <- ElPiGraph.R::computeElasticPrincipalGraph(
    Data = ScaffoldTree$CellCoordinates, NumNodes = N_yk_limits[1],
    NodesPositions = TopologyCoords, Edges = TopologyEdges,
    Method = 'CurveConfiguration', EP=lambda, RP=mu)

  for( i in N_yk_limits[2:length(N_yk_limits)])
  {
    InitCoords=ElasticTree[[1]]$Nodes
    InitEdges=ElasticTree[[1]]$Edges

    mu=(i-1)*mu_0
    lambda=((i-2)**3)*lambda_0

    ElasticTree <- ElPiGraph.R::computeElasticPrincipalGraph(
      Data = ScaffoldTree$CellCoordinates, NumNodes = i,
      NodesPositions = InitCoords, Edges = InitEdges,
      Method = 'CurveConfiguration', EP=lambda, RP=mu)
  }

  # Unlist the ElasticTree structure
  ElasticTree=ElasticTree[[1]]

  # Add the topology in the tree structure
  ETreeTopology= list(Endpoints=seq(from=1, to=length(ScaffoldTree$Endpoints), by=1), Branchpoints=seq(from=(length(ScaffoldTree$Endpoints) +1), to= (length(ScaffoldTree$Endpoints)) + length(ScaffoldTree$Branchpoints), by=1))
  ElasticTree=c(ElasticTree, Topology=1)
  ElasticTree$Topology <- ETreeTopology

  # Coords for the topology elements
  ElasticTree=c(ElasticTree, EndpointsCoords=1)
  ElasticTree$EndpointsCoords=ElasticTree$Nodes[ElasticTree$Topology$Endpoints,]
  ElasticTree=c(ElasticTree, BranchpointsCoords=1)
  ElasticTree$BranchpointsCoords=ElasticTree$Nodes[ElasticTree$Topology$Branchpoints,]

  # Add Input Coordinates to the tree structure
  ElasticTree=c(ElasticTree, CellCoords=1)
  ElasticTree$CellCoords <- ScaffoldTree$CellCoordinates

  # Add Tree Connectivity
  ElasticTree=c(ElasticTree, Connectivity=1)
  ElasticTree$Connectivity  <- Edges

  # Assign nodes from the embedded tree to the different branches and save the structure in the tree structure
  BranchesNodes=list()

  if(length(ScaffoldTree$Endpoints)==2)
  {
    path=TopologyEdges[1,1]
    first=TopologyEdges[1,1]
    while (first!=TopologyEdges[1,2])
    {
      edge=ElasticTree$Edges[which(ElasticTree$Edges[,1]==first),]
      second=edge[2]
      path=c(path, second)
      first=second
    }
    BranchesNodes[[1]] <- path
  }else
  {
    for (i in 1: dim(TopologyEdges)[1])
    {
      path=TopologyEdges[i,1]
      first=TopologyEdges[i,1]
      while (first!=TopologyEdges[i,2])
      {
        edge=ElasticTree$Edges[which(ElasticTree$Edges[,1]==first),]
        second=edge[2]
        path=c(path, second)
        first=second
      }
      BranchesNodes[[i]] <- path
    }
  }

  cell2yk=c()

  if(length(ScaffoldTree$Endpoints)==2)
  {
    cell_i=matrix(ScaffoldTree$CellCoordinates[1,], nrow=1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, ElasticTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk=rbind(cell2yk, c(1, closest_yk))
  }else
  {
    for (i in 1:dim(ScaffoldTree$CellCoordinates)[1])
    {
      cell_i=matrix(ScaffoldTree$CellCoordinates[i,], nrow=1)
      dist_cell_i=as.matrix(stats::dist(rbind(cell_i, ElasticTree$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
      #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
      closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
      cell2yk=rbind(cell2yk, c(i, closest_yk))
    }
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

#' Duplicate Elastic Tree Nodes
#'
#'Introduces intermediate nodes in between nodes being part of an edge in the ElasticTree structure
#' @param ElasticTree Elastic Tree to which nodes will be added
#' @return ElasticTree
#' @export
#'
DuplicateTreeNodes <- function(ElasticTree)
{
  ElasticTree2=ElasticTree
  N_yk=dim(ElasticTree$Nodes)[1]
  NewEdges=c()
  NewNodes=c()

  for(i in 1:dim(ElasticTree$Edges)[1])
    # for(i in 1:3)
  {
    Edge_i=ElasticTree$Edges[i,]
    coords_node1=ElasticTree$Nodes[Edge_i[1],]
    coords_node2=ElasticTree$Nodes[Edge_i[2],]
    coords_interpolate=coords_node1 + ((coords_node2 - coords_node1)/2)

    # node N_yk + i
    NewEdges=rbind(NewEdges, c(Edge_i[1], N_yk+i))
    NewEdges=rbind(NewEdges, c(Edge_i[2], N_yk+i))
    NewNodes=rbind(NewNodes, coords_interpolate)


    # Find branch
    branch_edges=c()
    for(j in 1:length(ElasticTree2$Branches))
    {
      if(Edge_i[1] %in% ElasticTree2$Branches[[j]] && Edge_i[2] %in% ElasticTree2$Branches[[j]])
      {
        branch_edges=j
      }
    }

    if(which(ElasticTree2$Branches[[branch_edges]]==Edge_i[2]) < which(ElasticTree2$Branches[[branch_edges]]==Edge_i[1]))
    {
      Edge_i=rev(Edge_i)
    }

    ElasticTree2$Branches[[branch_edges]]=c(ElasticTree2$Branches[[branch_edges]][1:which(ElasticTree2$Branches[[branch_edges]]==Edge_i[1])], N_yk+i, ElasticTree2$Branches[[branch_edges]][which(ElasticTree2$Branches[[branch_edges]]==Edge_i[2]):length(ElasticTree2$Branches[[branch_edges]])])
    # # add node to branch structure
    # if(!Edge_i[1] %in% ElasticTree$Topology$Branchpoints)
    # {
    #   j=1
    #   while(!Edge_i[1] %in% ElasticTree$Branches[[j]])
    #   {
    #     j=j+1
    #   }
    #   ElasticTree2$Branches[[j]]=c(ElasticTree2$Branches[[j]][1:which(ElasticTree2$Branches[[j]]==Edge_i[1])], N_yk+i, ElasticTree2$Branches[[j]][which(ElasticTree2$Branches[[j]]==Edge_i[2]):length(ElasticTree2$Branches[[j]])])
    # }else{
    #   j=1
    #   while(!Edge_i[2] %in% ElasticTree$Branches[[j]])
    #   {
    #     j=j+1
    #   }
    #   ElasticTree2$Branches[[j]]=c(ElasticTree2$Branches[[j]][1:which(ElasticTree2$Branches[[j]]==Edge_i[1])], N_yk+i, ElasticTree2$Branches[[j]][which(ElasticTree2$Branches[[j]]==Edge_i[2]):length(ElasticTree2$Branches[[j]])])
    #   # ElasticTree2$Branches[[j]]=c(ElasticTree2$Branches[[j]], N_yk+i)
    # }
  }

  ElasticTree2$Edges=NewEdges
  rownames(NewNodes)=c()
  ElasticTree2$Nodes=rbind(ElasticTree2$Nodes, NewNodes)

  # Reassing cells to nodes
  # reassigning cells to nodes on the full dimensional space
  cell2yk_post=c()
  for (i in 1:dim(ElasticTree2$CellCoords)[1])
  {
    cell_i= matrix(ElasticTree2$CellCoords[i,], nrow=1)
    dist_cell_i=as.matrix(stats::dist(rbind(cell_i, ElasticTree2$Nodes), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    closest_yk=sort(dist_cell_i[,1], index.return=T)$ix[2]-1
    cell2yk_post=rbind(cell2yk_post, c(i, closest_yk))
  }

  # allnodes=1:N_yk
  # cells_branchs_assigments=c()
  # for (i in 1:length(ElasticTree2$Branches))
  # {
  #   cells_branchs_assigments[which(cell2yk_post[,2] %in% ElasticTree2$Branches[[i]])]=i
  # }
  ElasticTree2$Cells2TreeNodes=cell2yk_post
  # --------------------

  return(ElasticTree2)
}

#' Get Cells from Trajectory
#'
#' Selects the cells that belong to a trajectory in an Elastric Tree between two elements in the tree topology.
#' @param ElasticTree Elastic tree from which cells will be taken
#' @param Start element in the topology from where the trajectory starts
#' @param End element in the topology where the trajectory ends
#' @param Pseudotimes pseudotime object to provide also the pseudotimes for the cells and nodes in the trajectory
#' @return Trajectory
#' @export
#' @importFrom igraph get.shortest.paths graph_from_adjacency_matrix
GetTrajectory <- function(ElasticTree, Start, End, Pseudotimes=NULL)
{
  # TrajectoryEPI=GetTrajectory(ElasticTree, Start=ElasticTree$Topology$Endpoints[1], End=ElasticTree$Topology$Endpoints[4], Pseudotimes = Pseudotimes)
  # Testing
  # Start=ElasticTree$Topology$Endpoints[1]
  # End=ElasticTree$Topology$Endpoints[2]
  # ElasticTree=ElasticTree
  # End Testing

  TreeTopologyNodes=c(ElasticTree$Topology$Endpoints, ElasticTree$Topology$Branchpoints)

  ElasticTree$Connectivity

  EdgesTopology=matrix(0, length(TreeTopologyNodes), length(TreeTopologyNodes))
  for(i in 1: dim(ElasticTree$Connectivity)[1])
  {
    EdgesTopology[ ElasticTree$Connectivity[i, 1],  ElasticTree$Connectivity[i, 2]]=1
    EdgesTopology[ ElasticTree$Connectivity[i, 2],  ElasticTree$Connectivity[i, 1]]=1
  }
  Graph_yk=igraph::graph_from_adjacency_matrix(EdgesTopology)

  path_brach_i=igraph::get.shortest.paths(Graph_yk,from = Start, to = End)
  TrajectoryBranchesPath=path_brach_i$vpath[[1]]

  NodesTrajectory=c()
  CellsTrajectory=c()

  # Check which branches are needed in the trajectory and order them in the right direction
  for(j in 1:(length(TrajectoryBranchesPath)-1))
  {
    print(paste(TrajectoryBranchesPath[j], TrajectoryBranchesPath[j+1]))
    for(i in 1:dim(ElasticTree$Connectivity)[1])
    {
      if(TrajectoryBranchesPath[j]==ElasticTree$Connectivity[i,1] && TrajectoryBranchesPath[j+1]==ElasticTree$Connectivity[i,2])
        NodesTrajectory=c(NodesTrajectory, ElasticTree$Branches[[i]])
      else if (TrajectoryBranchesPath[j]==ElasticTree$Connectivity[i,2] && TrajectoryBranchesPath[j+1]==ElasticTree$Connectivity[i,1])
        NodesTrajectory=c(NodesTrajectory, rev(ElasticTree$Branches[[i]]))
    }
  }

  # delete repeated nodes because of the branch pasting
  NodesTrajectory=unique(NodesTrajectory)

  NodesTrajectoryPseudotimes=1:length(NodesTrajectory)
  CellsTrajectoryPseudotimes=c()

  for(i in 1:length(NodesTrajectory))
  {
    CellsTrajectory=c(CellsTrajectory, ElasticTree$Cells2TreeNodes[which(ElasticTree$Cells2TreeNodes[,2]==NodesTrajectory[i]),1])
    CellsTrajectoryPseudotimes=c(CellsTrajectoryPseudotimes, rep(i, length(which(ElasticTree$Cells2TreeNodes[,2]==NodesTrajectory[i]))))
  }

  Trajectory= list(NodesTrajectory=NodesTrajectory, CellsTrajectory=CellsTrajectory, NodesTrajectoryPseudotimes=NodesTrajectoryPseudotimes, CellsTrajectoryPseudotimes=CellsTrajectoryPseudotimes)

  return(Trajectory)
}
