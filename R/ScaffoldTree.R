#' Calculate Scaffold tree
#'
#' creates the scaffold tree from the expression matrix in the Dataset object.
#' @param CellCoordinates file containing coordinates from the low dimensional manifold (e.g: first 3 diffusion components from a diffusion map)
#' @param NEndpoints Users can specify how many endpoints they want the algorithm to find. In case this variable is not defined all branches producing branches longer than sqrt(N/2) will be added to the tree structure
#' @param BranchMinLength Minimum number of nodes a branch has to contain in order to be included in the 3 structure. By default this value is set to sqrt(N/2) with N being the total number of cells in the dataset.
#' @param python_location url to the python3 executable binary. In case it is not specified a default call to python3 will be used.
#' @return ScaffoldTre object with the structure and connectivity of the Scaffold Tree
#' @export

CalculateScaffoldTree <- function(CellCoordinates, NEndpoints=NULL, BranchMinLength=-1, BranchMinLengthSensitive=-1, python_location="python3")
{
  CellCoordinates=as.matrix(CellCoordinates)
  CoordinatesFile=tempfile()
  write.table(CellCoordinates, file = CoordinatesFile, sep="\t", col.names = F, row.names = F)
  BranchMinLengthSensitive=floor(BranchMinLengthSensitive)
  ScaffoldTreeScript=paste(find.package("merlot"), "/python/ScaffoldTree.py", sep="")

  python_location <- glue::glue(
    "cd {find.package('merlot')}/venv",
    "source bin/activate",
    python_location,
    .sep = ";"
  )

  if(BranchMinLengthSensitive==-1 && is.null(NEndpoints))
  {
    BranchMinLengthSensitive=round(sqrt(dim(CellCoordinates)[1]))
  }

  if(is.null(NEndpoints))
  {
    #-------------------------------------------Execute TreeTopology.py-------------
    commands <- paste(python_location, " ",ScaffoldTreeScript, CoordinatesFile, "-BranchMinLength ", BranchMinLength, "-BranchMinLengthSensitive", BranchMinLengthSensitive)
  }  else
  {
    #-------------------------------------------Execute TreeTopology.py-------------
    commands <- paste(python_location, " ", ScaffoldTreeScript, CoordinatesFile, " -NBranches ", NEndpoints)
  }

  system(commands)

  # --------Read the topology elements from the TreeTopology.py output---------------
  ScaffoldTree=read_topology(CoordinatesFile, CellCoordinates)

  return(ScaffoldTree)
}

read_topology <-function (DataFile, CellCoordinates)
{
  TopologyData=read.table(file=paste(DataFile, "_TreeTopology.dat", sep=""), sep="\t", header=F, stringsAsFactors = F)
  # ----- We add 1 to the vectors because python numbers indexes from 0 instead of 1

  Endpoints=as.integer(unlist(strsplit(x=TopologyData[1,3], split=" "))) + 1

  if(length(Endpoints)==2)
  {
    Branchpoints=0
  }  else
  {
      Branchpoints=as.integer(unlist(strsplit(x=TopologyData[2,3], split=" "))) + 1
  }


  DijkstraPredecesors=read.table(file=paste(DataFile, "_DijkstraPredecesors.dat", sep=""), sep=" ", header=F, stringsAsFactors = F)
  DijkstraPredecesors= DijkstraPredecesors + 1
  DijkstraDistances=read.table(file=paste(DataFile, "_DijkstraDistances.dat", sep=""), sep=" ", header=F, stringsAsFactors = F)
  DijkstraSteps=read.table(file=paste(DataFile, "_DijkstraSteps.dat", sep=""), sep=" ", header=F, stringsAsFactors = F)

  # Reading the branches information
  Branches=as.matrix(TopologyData[3:dim(TopologyData)[1], 2:3])
  Branches=apply(Branches, 2, as.numeric)
  Branches=Branches+1
  rownames(Branches)=c()

  # Nodes that constitute the skeleton for the tree are saved in a vector
  SkeletonNodes=c()
  SkeletonEdges=c()
  BranchesNodes=list()
  BranchesCells=list()

  if(length(Endpoints)==2)
  {
    iBranch=calculate_path(Branches[1], Branches[2], DijkstraPredecesors)
    for(j in 1:(length(iBranch)-1))
    {
      SkeletonEdges=rbind(SkeletonEdges, (c(iBranch[j], iBranch[j+1])))
    }

    # Add nodes to the skeleton
    SkeletonNodes=c(SkeletonNodes, iBranch)
    # Add branch to the branch object
    BranchesNodes[[1]] <- iBranch
    names(BranchesNodes)=c()
    # Add all cells to the only branch
    BranchesCells[[1]] <- seq(1, dim(CellCoordinates)[1], 1)

    SkeletonNodes=sort(unique(SkeletonNodes))

    # Joining all the elements together in a list
    ScaffoldTree= list(Endpoints=Endpoints, Branchpoints=Branchpoints, DijkstraPredecesors=DijkstraPredecesors, DijkstraSteps=DijkstraSteps, DijkstraDistances= DijkstraDistances, Branches= Branches, SkeletonNodes=SkeletonNodes, SkeletonEdges= SkeletonEdges, CellCoordinates=CellCoordinates)
  }
  else if(length(Endpoints)>2)
  {
    for(i in 1:dim(Branches)[1])
    {
      iBranch=calculate_path(Branches[i,1], Branches[i,2], DijkstraPredecesors)
      for(j in 1:(length(iBranch)-1))
      {
        SkeletonEdges=rbind(SkeletonEdges, (c(iBranch[j], iBranch[j+1])))
      }
      # Add nodes to the skeleton
      SkeletonNodes=c(SkeletonNodes, iBranch)
      # Add branch to the branch object
      BranchesNodes[[i]] <- iBranch
      names(BranchesNodes)=c()
    }

    # Create the branchescells object as a temporal copy of the branchesnodes
    BranchesCells=BranchesNodes

    # Assign cells to branches
    Cells2ScaffoldCells=c()
    ScaffoldCells=sort(unique(unlist(BranchesNodes)))
    for(i in 1:dim(CellCoordinates)[1])
    {
      cell_i=matrix(CellCoordinates[i,], nrow=1)
      # Calculate which is the closest cell in the scaffold to cell i
      dist_cell_i=as.matrix(dist(rbind(cell_i, CellCoordinates[ScaffoldCells,]), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
      #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
      min_cell_dist=ScaffoldCells[sort(dist_cell_i[,1], index.return=T)$ix[2]-1]
      min_branch_dist=sort(dist_cell_i[,1], index.return=T)$x[2]
      Cells2ScaffoldCells=rbind(Cells2ScaffoldCells, c(i, min_cell_dist))
    }
    colnames(Cells2ScaffoldCells)=c("cell", "scaffoldcell")

    for(i in 1:length(BranchesNodes))
    {
      BranchesCells[[i]]= which(Cells2ScaffoldCells[,2] %in% BranchesNodes[[i]])
    }


    # cells_to_assign=which(!(1:dim(CellCoordinates)[1] %in% unlist(BranchesNodes)))
    #
    # for( i in cells_to_assign)
    # {
    #   # take cell i
    #   cell_i=matrix(CellCoordinates[i,], nrow=1)
    #
    #   # calculate minimum distance to each branch
    #   min_branch_dist=c()
    #   min_cell_dist=c()
    #
    #   for(j in 1:length(BranchesNodes))
    #   {
    #     dist_cell_i=as.matrix(dist(rbind(cell_i, CellCoordinates[BranchesNodes[[j]],]), method = "euclidean", diag = FALSE, upper = TRUE, p = 2))
    #     #find the closest yk index. Decrease the index in 1, since the 1 element is the element itself
    #     min_cell_dist=c(min_cell_dist, BranchesNodes[[j]][sort(dist_cell_i[,1], index.return=T)$ix[2]-1])
    #     min_branch_dist=c(min_branch_dist, sort(dist_cell_i[,1], index.return=T)$x[2])[1]
    #   }
    #
    #   # to which branch cell i should be assigned?
    #   branch_assignment=which(min_branch_dist==min(min_branch_dist))
    #   BranchesCells[[branch_assignment]]=c(BranchesCells[[branch_assignment]], i)
    # }

    Cells2Branches=rep(0, dim(CellCoordinates)[1])
    for(i in 1:length(BranchesCells))
    {
      Cells2Branches[BranchesCells[[i]]]=i
    }

    SkeletonNodes=sort(unique(SkeletonNodes))

    # Joining all the elements together in a list
    ScaffoldTree= list(Endpoints=Endpoints, Branchpoints=Branchpoints, DijkstraPredecesors=DijkstraPredecesors, DijkstraSteps=DijkstraSteps, DijkstraDistances= DijkstraDistances, Branches= Branches, SkeletonNodes=SkeletonNodes, SkeletonEdges= SkeletonEdges, CellCoordinates=CellCoordinates, Nodes2BranchesAssignments=BranchesNodes, Cells2BranchesAssignments=BranchesCells, Cells2Branches=Cells2Branches)
  }

  return(ScaffoldTree)
}

# Reconstructs cell path between two cells using the DijsktraPredecesors matrix calculated by TreeTopology.py
calculate_path <- function(cell_i, cell_j, DijkstraPredecesors)
{
  path = c()
  k=cell_j
  while ((k != cell_i) && (k >= 0))
  {
    path=c(path, k)
    k = DijkstraPredecesors[cell_i, k]
  }
  path=c(path, cell_i)
  return(path)
}

