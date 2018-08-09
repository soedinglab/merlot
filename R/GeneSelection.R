#' @importFrom stats var
LogVMR <- function(x) {
  return(log(x = stats::var(x = exp(x = x) - 1) / mean(x = exp(x = x) - 1)))
}

ExpMean <- function(x) {
  return(log(x = mean(x = exp(x = x) - 1) + 1))
}

#' seurat normalize
#'
#' Normalizes an expression matrix following the procedure as implemented in the Seurat's package
#' @param ExprMatrix cells x genes expression matrix to be normalized
#' @param scale.factor scaling factor to be used for normalizing the data
#' @export

seurat_normalize <- function(ExprMatrix, scale.factor = 1e4) {
  row.sums <- rowSums(ExprMatrix)
  by_cell_X <- (ExprMatrix / row.sums) * scale.factor
  logMatrix=log(by_cell_X + 1)
  return(logMatrix)
}

#' Seurat variable genes, according to the Seurat's package
#'
#' Identifies genes that are outliers on a 'mean variability plot'. First, calculates
#' average expression and dispersion for each gene. Next, divides genes into
#' num.bin (deafult 20) bins based on
#' their average expression, and calculates z-scores for dispersion within each
#' bin. The purpose of this is to identify variable genes while controlling for
#' the strong relationship between variability and average expression.
#'
#' Exact parameter settings may vary empirically from dataset to dataset, and
#' based on visual inspection of the plot.
#' Setting the y.cutoff parameter to 2 identifies genes that are more than two standard
#' deviations away from the average dispersion within a bin. The default X-axis function
#' is the mean expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
#' calculations are not performed in log-space, but the results are reported in log-space -
#' see relevant functions for exact details.
#'
#' @param ExprMatrix cells x genes expression matrix to be normalized
#' @param bin.size The bin size to use for the gene histogram
#' @param x.low.cutoff Bottom cutoff on x-axis for identifying variable genes
#' @param x.high.cutoff Top cutoff on x-axis for identifying variable genes
#' @param y.cutoff Bottom cutoff on y-axis for identifying variable genes
#' @param y.high.cutoff Top cutoff on y-axis for identifying variable genes
#' @param num.bin Total number of bins to use in the scaled analysis (default is 20)
#'
#' @return ScaffoldTre object with the structure and connectivity of the Scaffold Tree
#' @export
#'
#' @importFrom Matrix Matrix
seurat_variable_genes <- function(ExprMatrix,
                                  bin.size = 1000,
                                  x.low.cutoff = 0.1,
                                  x.high.cutoff = 8,
                                  y.cutoff = 1,
                                  y.high.cutoff = Inf,
                                  num.bin = 20)
  {
  genes.use <- rownames(x = t(ExprMatrix))
  gene.mean <- rep(x = 0, length(x = genes.use))
  names(x = gene.mean) <- genes.use
  gene.dispersion <- gene.mean
  gene.dispersion.scaled <- gene.mean
  max.bin <- floor(x = length(x = genes.use) / bin.size) + 1

  data <- Matrix::Matrix(t(ExprMatrix), sparse=TRUE)
  for (i in 1:max.bin) {
    my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1)) + 1
    my.inds <- my.inds[my.inds <= length(x = genes.use)]
    genes.iter <- genes.use[my.inds]
    data.iter <- data[genes.iter, , drop = F]
    gene.mean[genes.iter] <- apply(X = data.iter, MARGIN = 1, FUN = ExpMean)
    gene.dispersion[genes.iter] <- apply(X = data.iter, MARGIN = 1, FUN = LogVMR)
  }
  gene.dispersion[is.na(x = gene.dispersion)] <- 0
  gene.mean[is.na(x = gene.mean)] <- 0
  data_x_bin <- cut(x = gene.mean, breaks = num.bin)

  mean_y <- tapply(X = gene.dispersion, INDEX = data_x_bin, FUN = mean)
  sd_y <- tapply(X = gene.dispersion, INDEX = data_x_bin, FUN = sd)
  gene.dispersion.scaled <- (gene.dispersion - mean_y[as.numeric(x = data_x_bin)]) /
    sd_y[as.numeric(x = data_x_bin)]
  gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0
  mv.df <- data.frame(gene.mean, gene.dispersion, gene.dispersion.scaled)

  pass.cutoff <- names(x = gene.mean)[which(
    x = (
      (gene.mean > x.low.cutoff) & (gene.mean < x.high.cutoff)
    ) &
      (gene.dispersion.scaled > y.cutoff) &
      (gene.dispersion.scaled < y.high.cutoff)
  )]

  return(pass.cutoff)
}
