% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DiffusionMap.R
\name{wrap.DiffusionMap}
\alias{wrap.DiffusionMap}
\title{From destiny: Wrapper to create a diffusion map of cells.}
\usage{
wrap.DiffusionMap(data = destiny:::stopifnot_distmatrix(distance),
  sigma = "local", k = destiny::find_dm_k(dataset_n_observations(data,
  distance) - 1L), n_eigs = min(20L, dataset_n_observations(data,
  distance) - 2L), density_norm = TRUE, ..., distance = c("euclidean",
  "cosine", "rankcor"), n_local = seq(to = min(k, 7L), length.out =
  min(k, 3L)), rotate = FALSE, censor_val = NULL,
  censor_range = NULL, missing_range = NULL, vars = NULL,
  verbose = !is.null(censor_range), suppress_dpt = FALSE)
}
\arguments{
\item{data}{Expression data to be analyzed and covariates. Provide \code{vars} to select specific columns other than the default: all double value columns. If \code{distance} is a distance matrix, \code{data} has to be a \code{\link{data.frame}} with covariates only.}

\item{sigma}{Diffusion scale parameter of the Gaussian kernel. One of \code{'local'}, \code{'global'}, a (\link[base]{numeric}) global sigma or a \link[destiny]{Sigmas} object.
When choosing \code{'global'}, a global sigma will be calculated using \code{\link[destiny]{find_sigmas}}. (Optional. default: \code{'local'})
A larger sigma might be necessary if the eigenvalues can not be found because of a singularity in the matrix}

\item{k}{Number of nearest neighbors to consider (default: a guess betweeen 100 and \eqn{n - 1}. See \code{\link[destiny]{find_dm_k}}).}

\item{n_eigs}{Number of eigenvectors/values to return (default: 20)}

\item{density_norm}{logical. If TRUE, use density normalisation}

\item{...}{Unused. All parameters to the right of the \code{...} have to be specified by name (e.g. \code{DiffusionMap(data, distance = 'cosine')})}

\item{distance}{Distance measurement method applied to \code{data} or a distance matrix/\code{\link[stats]{dist}}. For the allowed values, see \code{\link[destiny]{find_knn}}}

\item{n_local}{If \code{sigma == 'local'}, the \code{n_local}th nearest neighbor(s) determine(s) the local sigma}

\item{rotate}{logical. If TRUE, rotate the eigenvalues to get a slimmer diffusion map}

\item{censor_val}{Value regarded as uncertain. Either a single value or one for every dimension (Optional, default: censor_val)}

\item{censor_range}{Uncertainity range for censoring (Optional, default: none). A length-2-vector of certainty range start and end. TODO: also allow \eqn{2\times G} matrix}

\item{missing_range}{Whole data range for missing value model. Has to be specified if NAs are in the data}

\item{vars}{Variables (columns) of the data to use. Specifying NULL will select all columns (default: All floating point value columns)}

\item{verbose}{Show a progressbar and other progress information (default: do it if censoring is enabled)}

\item{suppress_dpt}{Specify TRUE to skip calculation of necessary (but spacious) information for \code{\link[destiny]{DPT}} in the returned object (default: FALSE)}
}
\value{
a DiffusionMap object
}
\description{
The provided data can be a double \link[base]{matrix} of expression data or a \link[base]{data.frame} with all non-integer (double) columns
being treated as expression data features (and the others ignored), an \link[Biobase:class.ExpressionSet]{ExpressionSet}, or a \link[SingleCellExperiment]{SingleCellExperiment}.
}
