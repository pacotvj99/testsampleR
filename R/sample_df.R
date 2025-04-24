#' Example data of a test set
#'
#' Fake data with where each row corresponds to some observation in the annotated
#' test set, for which we know its stratum, sampling probability and true label.
#' This dataset is meant to be used for examples and illustrations of the functions
#' in the package.
#'
#' @format A data frame with 1394 rows and 5 variables:
#' \describe{
#'   \item{id}{Identifier of each sampled observation}
#'   \item{score}{Predicted probaility of exhibiting the outcome}
#'   \item{strata}{Stratum of the sampled observation (based on score)}
#'   \item{Prob}{Sampling probability}
#'   \item{truth}{Annotated labels: whether the observation actually exhibits the outcome}
#' }
#'
#' @docType data
#'
#' @usage data(sample_df)
#'
#' @keywords datasets
#'
#' @examples
#' data(sample_df)
"sample_df"
