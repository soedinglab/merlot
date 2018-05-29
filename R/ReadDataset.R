#' ReadDataset
#'
#' Read the dataset from a file containing gene the expression matrix
#' @param file filename including the path to where the expression matrix is contained. The file needs to have the cells as rows and the genes as columns.
#' @param featuresfile optional file with features for the different cells in the dataset
#' @return Dataset an object with the Expression Matrix, cell names and gene names
#' @usage ReadDataset(file)
#' @export
#'
#' @importFrom utils read.table
ReadDataset <- function(file, featuresfile=c())
{
  X=utils::read.table(file=file, sep="\t", header = T, stringsAsFactors = F, check.names = F, row.names = 1)
  ExpressionMatrix=as.matrix(X)
  dim(ExpressionMatrix)
  Descriptions=rownames(ExpressionMatrix)
  GeneNames=colnames(ExpressionMatrix)
  # Read Features for coloring Data

  Dataset= list(ExpressionMatrix=ExpressionMatrix, Descriptions=Descriptions, GeneNames=GeneNames)

  return (Dataset)
}

#' LoadExprMatrix
#'
#' @param ExprMatrix matrix object with n rows as cells and m columns as genes.
#' @return Dataset an object with the Expression Matrix, cell names and gene names
#' @usage LoadExprMatrix(file)
#' @export
LoadExprMatrix <-function(ExprMatrix)
{
  ExpressionMatrix=as.matrix(ExprMatrix)
  Descriptions=rownames(ExpressionMatrix)
  GeneNames=colnames(ExpressionMatrix)

  Dataset= list(ExpressionMatrix=ExpressionMatrix, Descriptions=Descriptions, GeneNames=GeneNames)
  return (Dataset)
}
