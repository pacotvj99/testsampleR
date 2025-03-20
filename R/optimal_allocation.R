#' Estimate optimal allocation
#'
#' This function shows how many observations to sample per stratum under optimal two-bin allocation
#'
#' This function can be used to determine how many observations to sample per stratum
#' under optimal allocation. It first determines how many positive observations to
#' sample in the test set, given parameter values and objective function, calling
#' \code{\link{optimal_n_positives}} internally. This is the value that minimizes an
#' objective function subject to a sample size constraint (see documentation of function
#' \code{\link{optimal_n_positives}} for details). Then, the function performs
#' proportional allocation within the positive and negative subsamples separately,
#' adjusting the number of sampled observations per bin so that the overall sample
#' size is as close to the target as possible.
#'
#' @param data A data.frame from which to sample observations. One column should
#' contain the strata used for sampling.
#' @param N_sample A numeric value capturing the desired sample size.
#' @param strata A character value capturing the column name of the sampling strata.
#' Column should contain factor or character values.
#' @param stratifying A character value capturing the column name of the stratifying
#' variable. This column should contain continuous numeric values between 0 and 1
#' (predicted probability of exhibiting the outcome).
#' @param min_per_bin A numeric value capturing the minimum number of observations per bin.
#' Default is 1.
#' @param threshold A numeric value capturing the threshold for binary classification.
#' @param pi1 A numeric value capturing the expected precision (share of cases
#' exhibiting the outcome out of predicted positives). This value must be entered
#' for the function to produce an output.
#' @param pi0 A numeric value capturing the inverse of the expected precision for
#' the negative class (share of cases exhibiting the outcome out of predicted
#' negatives). If left blank, the function computes `pi0` from recall and imbalance
#' if `pi0` is not explicitly entered.
#' @param recall A numeric value capturing the expected value of recall (share of
#' positive cases among those exhibiting the outcome). This value only needs to be
#' entered if `pi0` is left blank.
#' @param external_k A numeric value capturing the imbalance ratio (number of predicted
#' positives over predicted negatives) in the population on which the recall estimate
#' is based. This only needs to be entered if `recall` is based on a population
#' that has a different imbalance than the test's population. The argument is
#' only required if `pi0` is not entered directly. If left blank, it will be
#' calculated internally using `external_positive_share`, if this is entered, or
#' assumed to be equal to the test population's imbalance.
#' @param external_positive_share A numeric value capturing the positive share
#' (number of predicted positives out of all observations) in the population on
#' which the recall estimate is based. This only needs to be entered if `recall`
#' is based on a population that has a different imbalance than the test's population.
#' The argument is only required if `pi0` is not entered directly. If left blank,
#' it will be calculated internally using `external_positive_share`, if this is
#' entered, or assumed to be equal to the test population's imbalance.
#' @param weight_se_f1 A numeric value between 0 and 1 capturing how much to weigh
#' the SE for the F1 score in the objective function to be minimized. The default
#' for `weight_se_f1` is `1`.
#' @param weight_se_rec A numeric value between 0 and 1 capturing how much to weigh
#' the SE for recall in the objective function to be minimized. The default for
#' `weight_se_rec` is `0`.
#' @param weight_se_prec A numeric value between 0 and 1 capturing how much
#' to weigh the SE for precision in the objective function to be minimized. The
#' default for `weight_se_prec` is `0`.
#' @return A named vector whose values indicate the number of observations to be
#' sampled per each stratum, and whose names correspond to the values in the
#' `strata` variable of the input dataset.
#' @examples
#' ## example 1: say that we want to sample observations using two-bin stratification
#' ## of positives and negatives (dichotomized at 0.5 threshold), picking the number
#' ## of positives and negatives optimally, subject to a sample size constraint that
#' ## N_positive + N_negative = 500. here, assume we're minimizing the objective
#' ## function SE(F1) + 0.5 SE(precision) + 0.5 SE(recall). as assumptions, we have
#' ## guesses for pi1 (TP/TP+FP) and pi0 (FN/FN+TN) of 0.4 and 0.02.
#' data(pop_df)
#' pop_df$strata <- ifelse(pop_df$score>=0.5, "yes", "no")
#' optimal_allocation(data = pop_df, N_sample = 500, strata = 'strata',
#'                    stratifying = 'score', min_per_bin=1, threshold=0.5, pi1=0.4,
#'                    pi0=0.2, weight_se_f1=1, weight_se_rec=0.5, weight_se_prec=0.5)
#' ## example 2: same but we now want to further implement proportional stratification
#' ## within the positives and negatives separately, after having sampled the optimal
#' ## number of positives and negatives subject to the sample size constraint. also,
#' ## instead of making a guess about pi0, we make a guess about recall (0.6) based on
#' ## previous literature, and specify the imbalance in the data were said recall was
#' ## estimated (20% positives). objective function and constraints are as above.
#' data(pop_df)
#' pop_df$strata <- cut(pop_df$score, breaks=seq(0,1,0.1))
#' optimal_allocation(data = pop_df, N_sample = 500, strata = 'strata',
#'                    stratifying = 'score', min_per_bin=1, threshold=0.5, pi1=0.4,
#'                    recall=0.6, external_positive_share=0.2,
#'                    weight_se_f1=1, weight_se_rec=0.5, weight_se_prec=0.5)
#' ## example 3: like example 2 but assuming only care about F1, and assuming we
#' ## specify the imbalance ratio in the dataset where recall was estimated: we
#' ## know it has 1 positive for every 2 negatives so external_k is 1/2=0.5
#' data(pop_df)
#' pop_df$strata <- cut(pop_df$score, breaks=seq(0,1,0.1))
#' optimal_allocation(data = pop_df, N_sample = 500, strata = 'strata',
#'                    stratifying = 'score', min_per_bin=1, threshold=0.5,
#'                    pi1=0.4, recall=0.6, external_k=0.5,
#'                    weight_se_f1=1, weight_se_rec=0, weight_se_prec=0)
#' ## note the function disregards NAs in the strata variable
#' ## note that for optimal allocation the strata must be based somehow on predicted
#' ## probabilities (quantile- or fixed-intervals, where each stratum contains only
#' ## positives or only negatives but not both)
#' @export
optimal_allocation <- function(data, N_sample, strata, stratifying, min_per_bin=1, threshold=0.5, pi1=NULL,
                               pi0=NULL, recall=NULL, external_k=NULL, external_positive_share=NULL,
                               weight_se_f1=1, weight_se_rec=0, weight_se_prec=0){
  if("data.frame" %in% class(data)){
    class(data) <- "data.frame"
  } else {
    stop("Invalid input data")
  }
  if(!is.character(strata)){
    stop("Invalid strata argument")
  }
  if(!is.character(stratifying)){
    stop("Invalid stratifying argument")
  }

  if(sum(colnames(data)==stratifying)!=1){
    stop("Column with stratifying variable cannot be uniquely identified")
  }
  if(sum(colnames(data)==strata)!=1){
    stop("Column with strata cannot be uniquely identified")
  }
  if(any(is.na(data$stratifying))){
    warning("NAs in stratifying variable")
  }
  if(any(is.na(data$strata))){
    warning("NAs in strata variable")
  }
  data$stratifying <- unlist(data[,stratifying])
  data$strata <- unlist(data[,strata])

  if(threshold>1 | threshold<0){
    stop("Invalid threshold")
  }
  if(!is.numeric(N_sample) || !is.numeric(min_per_bin)){
    stop("Invalid binning specifications")
  } else if(((N_sample%%1)!=0) || ((min_per_bin%%1)!=0)){
    stop("Invalid binning specifications")
  }

  if(is.factor(data$strata)){
    data$strata <- as.character(data$strata)
  } else if(!is.character(data$strata)){
    stop("Invalid strata variable")
  }
  if(!is.numeric(data$stratifying)){
    stop("Stratifying variable is not numeric")
  } else if(min(data$stratifying)<0 | max(data$stratifying)>1){
    stop("Stratifying variable is numeric but not probability")
  } else if(all(data$stratifying>=threshold) | all(data$stratifying<threshold)){
    stop("All stratifying values fall to same side of the threshold")
  }

  #stratify sample into positives and negatives
  sam_pos <- data[data$stratifying>=threshold,]
  sam_neg <- data[data$stratifying<threshold,]
  if(any(sam_pos$strata %in% unique(sam_neg$strata)) |
     any(sam_neg$strata %in% unique(sam_pos$strata))){
    stop("All stratifying values fall to same side of the threshold")
  }

  #estimate optimal number of positives
  k <- sum(data$stratifying>=threshold)/sum(data$stratifying<threshold)
  optimal_n <- optimal_n_positives(N_sample, pi1, pi0, recall, k=k, positive_share=NULL,
                                   external_k, external_positive_share,
                                   weight_se_f1, weight_se_rec, weight_se_prec)
  if((!is.numeric(optimal_n))|(length(optimal_n)>1)){
    stop("Invalid optimal number of positives to be sampled")
  }

  #if optimally sample more positives than there are, coerce the sampled number of positives to the number of positives. same for negatives
  if(nrow(sam_pos)<optimal_n){
    optimal_n <- nrow(sam_pos)
    warning("Data contains fewer positives than optimally sampled")
  } else if(nrow(sam_neg)<(N_sample - optimal_n)){
    optimal_n <- N_sample - nrow(sam_neg)
    warning("Data contains fewer negatives than optimally sampled")
  }

  #optimal allocation for positives and negatives separately
  myallocation_neg <- proportional_allocation(sam_neg, (N_sample - optimal_n), strata, min_per_bin)
  myallocation_pos <- proportional_allocation(sam_pos, optimal_n, strata, min_per_bin)
  myallocation <- c(myallocation_neg, myallocation_pos)
  return(myallocation)
}


#' Calculator of optimal number of positives
#'
#' This function estimates how many positives should optimally be sampled
#'
#' This function can be used to estimate how many positive observations should be
#' included in the test set under efficient stratified random sampling (SRS). This
#' is the value that minimizes an objective function subject to a sample size constraint.
#' The objective function is a linear combination of the standard errors (SEs) of
#' precision, recall and F1 score, where weights can be specified by the user.
#' These standard errors (SEs) are based on analytic formulas for the
#' variance of the sampling distribution of the test set estimates of F1, precision
#' and recall under two-bin stratified sampling. Such variances are based on expected
#' parameter values, which need to be entered by the user. The value returns a
#' unique integer value that minimizes the expression subject to the constraint.
#'
#' @param N_sample A numeric value capturing the total test set sample size considered
#' as binding constraint (includes both positives and negatives).
#' @param pi1 A numeric value capturing the expected precision (share of cases
#' exhibiting the outcome out of predicted positives). This value must be entered
#' for the function to produce an output.
#' @param pi0 A numeric value capturing the inverse of the expected precision for
#' the negative class (share of cases exhibiting the outcome out of predicted
#' negatives). If left blank, the function computes `pi0` from recall and imbalance
#' if `pi0` is not explicitly entered.
#' @param recall A numeric value capturing the expected value of recall (share of
#' positive cases among those exhibiting the outcome). This value only needs to be
#' entered if `pi0` is left blank.
#' @param k A numeric value capturing the imbalance ratio (number of predicted positives
#' over predicted negatives) in the population from which the test set should be
#' drawn. Argument `k` can be left blank, but then `positive_share` must be specified.
#' @param positive_share A numeric value capturing the positive share (number of
#' predicted positives out of all observations) in the population from which the
#' test set should be drawn. Argument `positive_share` can be left blank, but
#' then `k` must specified.
#' @param external_k A numeric value capturing the imbalance ratio (number of predicted
#' positives over predicted negatives) in the population on which the recall estimate
#' is based. This only needs to be entered if `recall` is based on a population
#' that has a different imbalance than the test's population. The argument is
#' only required if `pi0` is not entered directly. If left blank, it will be
#' calculated internally using `external_positive_share`, if this is entered, or
#' assumed to be equal to the test population's imbalance.
#' @param external_positive_share A numeric value capturing the positive share
#' (number of predicted positives out of all observations) in the population on
#' which the recall estimate is based. This only needs to be entered if `recall`
#' is based on a population that has a different imbalance than the test's population.
#' The argument is only required if `pi0` is not entered directly. If left blank,
#' it will be calculated internally using `external_positive_share`, if this is
#' entered, or assumed to be equal to the test population's imbalance.
#' @param weight_se_f1 A numeric value between 0 and 1 capturing how much to weigh
#' the SE for the F1 score in the objective function to be minimized by the calculator.
#' The default for `weight_se_f1` is `1`.
#' @param weight_se_rec A numeric value between 0 and 1 capturing how much to weigh
#' the SE for recall in the objective function to be minimized by the calculator.
#' The default for `weight_se_rec` is `0`.
#' @param weight_se_prec A numeric value between 0 and 1 capturing how much
#' to weigh the SE for precision in the objective function to be minimized by the
#' calculator. The default for `weight_se_prec` is `0`.
#' @return A numeric value capturing how many positive observations should be sampled
#' given the parameter values and the objective function specified by the user.
#' @export
optimal_n_positives <- function(N_sample=NULL, pi1=NULL, pi0=NULL, recall=NULL, k=NULL,
                                positive_share=NULL, external_k=NULL, external_positive_share=NULL,
                                weight_se_f1=1, weight_se_rec=0, weight_se_prec=0){
  if(!is.numeric(k) & is.numeric(positive_share)){
    k <- get_k(positive_share)
  } else if(is.numeric(k) & is.numeric(positive_share)){
    message("Both k and positive share specified: k used for analysis")
  }
  if(!is.numeric(external_k) & is.numeric(external_positive_share)){
    external_k <- get_k(external_positive_share)
  } else if(is.numeric(external_k) & is.numeric(external_positive_share)){
    message("Both external k and external positive share specified: external k used for analysis")
  }
  if(is.numeric(pi0) & is.numeric(pi1) & is.numeric(k)){
    implied_rec <- 1/(1+(pi0/(k*pi1)))
    if(implied_rec<0.02){
      message("Very small recall implied for test set: consider revising assumptions")
    }
  }
  if(!is.numeric(external_k) & !is.numeric(pi0) & is.numeric(recall)){
    message("In-sample imbalance assumed to estimate pi0")
    external_k <- k
  } else if(!is.numeric(pi0) & !is.numeric(recall)){
    message("Both pi0 and recall specified: only pi0 used for analysis")
  }
  if(!is.numeric(pi0) & is.numeric(pi1) & is.numeric(recall) & is.numeric(external_k)){
    pi0 <- get_pi0(pi1, recall, external_k)
  }
  if(!is.numeric(weight_se_f1) | !is.numeric(weight_se_rec) | !is.numeric(weight_se_prec) |
     !is.numeric(pi1) | !is.numeric(pi0) | !is.numeric(k) | !is.numeric(N_sample)){
    stop("Some parameters have not been rightly specified")
  }
  if(pi1<=0 | pi1>=1 | pi0<=0 | pi0>=1 | ((N_sample %% 1) !=0)){
    stop("Some parameters have not been rightly specified")
  }

  #estimate optimal number of positives
  recall <- k*pi1/(k*pi1+pi0)
  f1 <- 2*(recall*pi1)/(recall+pi1)
  f_star <- f1/(2-f1)

  f1_variance <- (4/(1+f_star)^4)*f_star^2*(((1-pi1)/(pi1*(1:(N_sample-1)))) + ((1-pi0)*pi0/((pi0+k)^2*(N_sample-(1:(N_sample-1))))))
  recall_variance <- (pi0/(k*pi1*((1+(pi0/(k*pi1)))^2)))^2*(((1-pi1)/(pi1*(1:(N_sample-1)))) + ((1-pi0)/(pi0*(N_sample-(1:(N_sample-1))))))
  precision_variance <- (pi1*(1-pi1)/((1:(N_sample-1))))
  optimal_n <- which.min(weight_se_f1*sqrt(f1_variance) + weight_se_rec*sqrt(recall_variance) + weight_se_prec*sqrt(precision_variance))

  return(optimal_n)
}
