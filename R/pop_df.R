#' Example data of a population from which to draw a test set
#'
#' Fake data where each row corresponds to some observation in the population from
#' which the test set is to be sampled, and for which we know its predicted
#' probability of exhibiting the outcome of interest according to some classifier.
#' This dataset is meant to be used for examples and illustrations of the functions
#' in the package.
#'
#' @format A data frame with 20000 rows and 2 variables:
#' \describe{
#'   \item{id}{Identifier of each sampled observation}
#'   \item{score}{Predicted probaility of exhibiting the outcome}
#' }
#'
#' @docType data
#'
#' @usage data(pop_df)
#'
#' @keywords datasets
#'
#' @examples
#' data(pop_df)
"pop_df"

