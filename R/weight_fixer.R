#' Fix weights after incomplete annotations
#'
#' This function estimates probability weights given a partially annotated test set
#'
#' This function can be used on a test set drawn by stratified sampling that has
#' not been fully annotated, e.g. if human annotators stopped annotating halfway
#' through the annotation task. The function re-computes sampling probabilities
#' given the annotation was incomplete. This assumes that the comments that were
#' indeed annotated from each stratum are a SRS of the stratum.
#'
#' @param data A data.frame containing the annotated test set.
#' @param truth A character value capturing the column name of the annotated labels.
#' Column should capture whether each observation exhibits the outcome (1/TRUE)
#' or not (0/FALSE). NAs should be entered for observations that were not annotated.
#' @param strata A character value capturing the column name of the sampling strata.
#' Column should contain factor or character values.
#' @param probs A character value capturing the column name of the sampling probabilities.
#' Default assumes equal sampling probability for all observations in the original sampling.
#' @return A data.frame just like the input `data`, but with recomputed sampling
#' probabilities in column `Prob_new`.
#' @examples
#' ## assume that some observations were not annotated, such that their value in
#' ## the column with annotated labels is NA. function will update the weights so
#' ## estimates remain unbiased. Note this assumes that the annotated observations
#' ## are a SRS of their respective stratum if the test set was drawn by stratified
#' ## sampling, or a SRS of the entire population if the test set was drawn by SRS.
#' data(sample_df)
#' sample_df$truth[30:45] <- NA #code some labels as missing because not annotated
#' weight_fixer(data = sample_df, truth = "truth", strata = "strata", probs="Prob")
#' @export
weight_fixer <- function(data, truth="truth", strata="strata", probs=NULL){
  if("data.frame" %in% class(data)){
    class(data) <- "data.frame"
  } else {
    stop("Invalid input data")
  }

  for(conflicting_column in c("Prob", "strata", "truth", "Prob_new")){
    if(!((truth=="truth" & conflicting_column=="truth") | (strata=="strata" & conflicting_column=="strata") | (probs=="Prob" & conflicting_column=="Prob"))){
      if(conflicting_column %in% colnames(data)){
        if(paste0(conflicting_column, ".0") %in% colnames(data)){
          stop(paste0("Invalid variable name in input dataset: ", conflicting_column, ".0"))
        } else {
          colnames(data)[colnames(data)==conflicting_column] <- paste0(conflicting_column, ".0")
          message(paste0("Column '", conflicting_column,"' renamed to '", paste0(conflicting_column, ".0"), "'"))
        }
      }
    }
  }

  if(is.null(probs)){
    warning("No probabilities in data: constant sampling probability assumed. This is a problem if stratified sampling was used for drawing test set")
    data$Prob <- 1
  } else if(!is.character(probs)){
    stop("Invalid probs argument")
  } else if(sum(colnames(data)==probs)!=1){
    stop("Column with probability labels cannot be uniquely identified")
  } else {
    data$Prob <- unlist(data[,probs])
  }
  if(any(data$Prob<0) || any(data$Prob>1)){
    stop("Invalid sampling probabilities")
  }

  if(!is.character(truth)){
    stop("Invalid truth argument")
  } else if(sum(colnames(data)==truth)!=1){
    stop("Column with annotated labels cannot be uniquely identified")
  } else {
    data$truth <- unlist(data[,truth])
  }

  if(!is.character(strata)){
    stop("Invalid strata argument")
  } else if(sum(colnames(data)==strata)!=1){
    stop("Column with stratum labels cannot be uniquely identified")
  } else {
    data$strata <- unlist(data[,strata])
  }


  if(any(is.na(data$truth))){
    # output <- data %>%
    #   dplyr::group_by(strata) %>%
    #   dplyr::mutate(number_of_sampled = dplyr::n(),
    #                 number_of_coded = sum(!is.na(truth))) %>%
    #   dplyr::ungroup() %>%
    #   dplyr::mutate(number_in_population=number_of_sampled/Prob) %>%
    #   dplyr::mutate(Prob_new=number_of_coded/number_in_population) %>%
    #   dplyr::pull(Prob_new)
    output <- within(data, {
      number_of_sampled <- as.numeric(as.character(ave(strata, strata, FUN = length)))
      number_of_coded <- as.numeric(ave(!is.na(truth), strata, FUN = sum))
      number_in_population <- number_of_sampled / Prob
      Prob_new <- number_of_coded / number_in_population
    })
    output <- output$Prob_new
  } else{
    message("No weight needs fixing")
    output <- data$Prob
  }
  return(output)
}
