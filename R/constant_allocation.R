#' Estimate constant allocation
#'
#' This function shows how many observations to sample per stratum under constant allocation
#'
#' This function can be used to determine how many observations to sample per stratum
#' under constant allocation. It first divides the desired sample size equally across
#' all strata, and then rounds the allocated observations per bin, and lastly adjusts
#' them so the overall sample size is as close to the target as possible.
#'
#' @param data A data.frame from which to sample observations. One column should
#' contain the strata used for sampling.
#' @param N_sample A numeric value capturing the desired sample size.
#' @param strata A character value capturing the column name of the sampling strata.
#' Column should contain factor or character values.
#' @param min_per_bin A numeric value capturing the minimum number of observations per bin.
#' Default is 1.
#' @return A named vector whose values indicate the number of observations to be
#' sampled per each stratum, and whose names correspond to the values in the
#' `strata` variable of the input dataset.
#' @examples
#' ## say that we want to sample a constant number of observations per stratum
#' ## but with appropriate rounding given one cannot sample fractional observations.
#' ## this function tells us how many obs per stratum to sample, given the strata
#' data(pop_df)
#' pop_df$strata <- cut(pop_df$score, 10)
#' constant_allocation(data=pop_df, N_sample=375, strata="strata", min_per_bin=1)
#' ## note the function disregards NAs in the strata variable
#' ## note the strata can be based on any information: either predicted probabilities
#' ## (quantile- or fixed-intervals) or any other categorical information
#' @export
constant_allocation <- function(data, N_sample, strata, min_per_bin=1){
  if("data.frame" %in% class(data)){
    class(data) <- "data.frame"
  } else {
    stop("Invalid input data")
  }
  if(!is.numeric(N_sample) || !is.numeric(min_per_bin)){
    stop("Invalid binning specifications")
  } else if(((N_sample%%1)!=0) || ((min_per_bin%%1)!=0)){
    stop("Invalid binning specifications")
  }

  if(!is.character(strata)){
    stop("Invalid strata argument")
  }
  if(sum(colnames(data)==strata)!=1){
    stop("Column with strata cannot be uniquely identified")
  }
  data$strata <- unlist(data[,strata])
  if(is.factor(data$strata)){
    data$strata <- as.character(data$strata)
  } else if(!is.character(data$strata)){
    stop("Invalid strata variable")
  }
  if(any(is.na(data$strata))){
    warning("NAs in strata variable")
  }

  #create count of observations per bin
  strata_count <- table(data$strata)
  strata_count <- strata_count[order(names(strata_count))]
  if((min_per_bin*length(strata_count))>N_sample){
    stop("Error: too many bins given the sample size")
  }
  strata_count <- strata_count[strata_count>0]

  #one first naive go at constant allocation
  myallocation <- rep(N_sample/length(strata_count), length(strata_count))
  names(myallocation) <- names(strata_count)

  any_issues <- TRUE
  iter <- 1
  while(any_issues && iter<1000){
    #if allocation implies sampling fewer than 1, coerce it to 1
    which_binding <- myallocation<min_per_bin
    myallocation <- ifelse(which_binding, min_per_bin, myallocation)
    #if allocation implies sampling more obs than there are, coerce it to the number of obs in the bin
    strata_count_off <- (strata_count - myallocation) <0
    myallocation[strata_count_off] <- strata_count[strata_count_off]
    #rebalance the allocation so it still sums to N_sample
    needup <- N_sample-sum(myallocation)
    myallocation[!which_binding & !strata_count_off] <- myallocation[!which_binding & !strata_count_off]+
      (needup/sum(!which_binding & !strata_count_off))
    #round up allocations so it is as close as possible to N_sample
    remainders <- myallocation %% 1
    if(any(remainders!=0)){
      myallocation <- round(myallocation)
      #rebalance after rounding so it still is as close as possible to N_sample
      needup <- N_sample-sum(myallocation)
      if(needup>0){
        myallocation[names(sort(remainders[remainders<=0.5], decreasing=TRUE))[1:needup]] <-
          myallocation[names(sort(remainders[remainders<=0.5], decreasing=TRUE))[1:needup]] + 1
      } else if(needup<0){
        myallocation[names(sort(remainders[remainders>=0.5]))[1:(-needup)]] <-
          myallocation[names(sort(remainders[remainders>=0.5]))[1:(-needup)]] - 1
      }
    }
    #if sample more obs than there are in some bin, if sample zero obs in some bin, or if some fewer than min_per_bin while possible to sample more, then keep going
    any_issues <- any((myallocation<min_per_bin & strata_count>=min_per_bin) | (myallocation<1) | (strata_count < myallocation))
    iter <- iter + 1
  }

  if(iter==1000){
    stop("Error: Valid allocation could not be found")
  }
  return(myallocation)
}
