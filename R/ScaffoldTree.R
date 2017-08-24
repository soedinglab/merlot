#' Calculate Scaffold tree
#'
#' creates the scaffold tree from the expression matrix in the Dataset object.
#' @param CellCoordinates file containing coordinates from the low dimensional manifold (e.g: first 3 diffusion components from a diffusion map)
#' @param NEndpoints Users can specify how many endpoints they want the algorithm to find. In case this variable is not defined all branches producing branches longer than sqrt(N/2) will be added to the tree structure
#' @param python_location url to the python3 executable binary. In case it is not specified a default call to python3 will be used.
#' @return ScaffoldTre object with the structure and connectivity of the Scaffold Tree
#' @export

CalculateScaffoldTree <- function(CellCoordinates, NEndpoints=NULL, python_location="python3")
{
  CoordinatesFile=tempfile()
  write.table(CellCoordinates, file = CoordinatesFile, sep="\t", col.names = F, row.names = F)

  ScaffoldTreeScript=paste(find.package("merlot"), "/python/ScaffoldTree.py", sep="")

  if(is.null(NEndpoints))
  {
    #-------------------------------------------Execute TreeTopology.py-------------
    system(paste(python_location, " ",ScaffoldTreeScript, CoordinatesFile), wait = TRUE)
  }  else
  {
    #-------------------------------------------Execute TreeTopology.py-------------
    system(paste(python_location, " ", ScaffoldTreeScript, CoordinatesFile, " -NBranches ", NEndpoints), wait = TRUE)
  }

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

  if(length(Endpoints)==2)
  {
    iBranch=calculate_path(Branches[1], Branches[2], DijkstraPredecesors)
    for(j in 1:(length(iBranch)-1))
    {
      SkeletonEdges=rbind(SkeletonEdges, (c(iBranch[j], iBranch[j+1])))
    }
    SkeletonNodes=c(SkeletonNodes, iBranch)
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
      SkeletonNodes=c(SkeletonNodes, iBranch)
    }

    SkeletonNodes=sort(unique(SkeletonNodes))

    # Joining all the elements together in a list
    ScaffoldTree= list(Endpoints=Endpoints, Branchpoints=Branchpoints, DijkstraPredecesors=DijkstraPredecesors, DijkstraSteps=DijkstraSteps, DijkstraDistances= DijkstraDistances, Branches= Branches, SkeletonNodes=SkeletonNodes, SkeletonEdges= SkeletonEdges, CellCoordinates=CellCoordinates)
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

