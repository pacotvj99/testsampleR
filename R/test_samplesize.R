#' Sample size calculator
#'
#' This function estimates how many observations need to be sampled to achieve a desired standard error.
#'
#' This function can be used to estimate how many observations need to be included
#' in the test set under simple random sampling (SRS) or under efficient stratified
#' sampling. The standard errors (SEs) are based on analytic formulas for the
#' variance of the sampling distribution of the test set estimates of F1, precision
#' and recall. Such variances are based on expected parameter values, which need to
#' be entered by the user. Given parameter values, the function looks at what sample
#' sizes up to `max_N` verify that the SEs for F1, precision and recall are below
#' the maximum desired SE, which need to be entered by the user. Under SRS, the
#' function returns the minimum sample size that verifies the desired SEs conditions,
#' as well as the SEs under such sample size. Under stratified sampling, the function
#' returns the minimum sample size such that the desired SE conditions can be verified
#' for some allocation (i.e. some number of positives to be sampled): this is based
#' on two-bin stratified sampling (sampling positives and negatives separately), and
#' the resulting sample calculation is not valid if stratifying on some other variable.
#' Under stratified sampling, the function also returns how to allocate the sample
#' size between positives and negatives to verify the desired SEs conditions for
#' the given sample size. When there are multiple values of `n_positive` that
#' verify the desired SE conditions for a given sample size, the value of `n_positive`
#' is selected by picking the number that minimizes the sum of the SE for F1 score,
#' precision and recall, weighted by some user-provided weights.
#'
#' @param se_f1 A numeric value capturing the desired standard error for the F1
#' score of the positive class. If non-numeric, then the sample size calculation
#' ignores the SE of F1. The sample size calculator selects a sample size such that
#' the SE of the F1 score is below `se_f1`.
#' @param se_precision A numeric value capturing the desired standard error for
#' precision of the positive class. If non-numeric, then the sample size calculation
#' ignores the SE of precision. The sample size calculator selects a sample size
#' such that the SE of precision is below `se_precision`.
#' @param se_recall A numeric value capturing the desired standard error for
#' recall of the positive class. If non-numeric, then the sample size calculation
#' ignores the SE of recall. The sample size calculator selects a sample size such
#' that the SE of recall is below `se_recall`.
#' @param max_N A numeric value capturing the maximum sample size to be considered.
#' @param min_N A numeric value capturing the minimum sample size to be considered.
#' Default is `1`.
#' @param by_N A numeric value capturing the size of the intervals between `min_N`
#' and `max_N` to be considered as possible sample sizes. Default is `1`.
#' @param pi1 A numeric value capturing the expected precision (share of cases
#' exhibiting the outcome out of predicted positives). This value must be entered
#' for the function to produce an output.
#' @param pi0 A numeric value capturing the inverse of the expected precision for
#' the negative class (share of cases exhibiting the outcome out of predicted
#' negatives). If left blank, the function computes `pi0` from recall and imbalance
#' if `pi0` is not explicitly entered. Note this value need not be entered if
#' `se_recall` and `se_f1` are left blank.
#' @param recall A numeric value capturing the expected value of recall (share of
#' positive cases among those exhibiting the outcome). This value only needs to be
#' entered if `se_recall` or `se_f1` are specified but `pi0` is left blank.
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
#' assumed to be equal to the test population's imbalance `k`.
#' @param external_positive_share A numeric value capturing the positive share
#' (number of predicted positives out of all observations) in the population on
#' which the recall estimate is based. This only needs to be entered if `recall`
#' is based on a population that has a different imbalance than the test's population.
#' The argument is only required if `pi0` is not entered directly. If left blank,
#' it will be calculated internally using `external_positive_share`, if this is
#' entered, or assumed to be equal to the test population's imbalance.
#' @param weight_se_f1 A numeric value between 0 and 1 capturing how much to weigh
#' the SE for the F1 score. This only applies when there are multiple values of
#' `n_positive` that verify the desired SE conditions for a given sample size under
#' stratified sampling. When this happens, the value of `n_positive` is selected by
#' picking the number that minimizes the sum of the SE for F1 score, precision and
#' recall, weighted by `weight_se_f1`, `weight_se_precision` and , `weight_se_recall`.
#' The default for `weight_se_f1` is `1`.
#' @param weight_se_prec A numeric value between 0 and 1 capturing how much to weigh
#' the SE for precision. This only applies when there are multiple values of `n_positive`
#' that verify the desired SE conditions for a given sample size under  stratified
#' sampling. When this happens, the value of `n_positive` is selected by picking
#' the number that minimizes the sum of the SE for F1 score, precision and recall,
#' weighted by `weight_se_f1`, `weight_se_prec` and , `weight_se_rec`. The default
#' for `weight_se_prec` is `1`.
#' @param weight_se_rec A numeric value between 0 and 1 capturing how much to
#' weigh the SE for recall. This only applies when there are multiple values of
#' `n_positive` that verify the desired SE conditions for a given sample size under
#' stratified sampling. When this happens, the value of `n_positive` is selected by
#' picking the number that minimizes the sum of the SE for F1 score, precision and
#' recall, weighted by `weight_se_f1`, `weight_se_prec` and , `weight_se_rec`. The
#' default for `weight_se_rec` is `1`.
#' @return A list of two data.frame objects (`srs` and `stratified`), which contain the result of the power
#' calculation under SRS and stratified sampling respectively. These data.frames
#' contain the minimum number of observations (`sample_size_srs` and `sample_size_strat`)
#' needed to verify the desired SEs conditions (i.e. smaller than `se_f1`, `se_precision`
#' and `se_recall`). For stratified sampling, the output also contains how many
#' positive observations (`n_positives`) need to be sampled to achieve the desired SEs.
#' Under either sampling, the output also indicates what the SE would be for each
#' metric and sampling stratgy, given the selected sample size. When the maximum
#' sample size (`max_N`) is too small to achieve the desired SEs, the output indicates
#' this instead.
#' @examples
#' ## example 1: required sample size to get SE of F1 below 3pp
#' ## assuming we have guesses for pi1 (TP/TP+FP) and pi0 (FN/FN+TN) of 0.4 and 0.02
#' ## assuming we know in the test set 10% of observations are positive
#' ## further assuming that when multiple allocations exist, we only care about F1
#' test_samplesize(se_f1=0.03, max_N=4000, pi1=0.4, pi0=0.02, positive_share=0.1,
#'                weight_se_f1=1, weight_se_prec=0, weight_se_rec=0)
#' ## example 2: same but assuming we know test set has 1 positive for every 2
#' ## negatives so k is 1/2=0.5
#' test_samplesize(se_f1=0.03, max_N=4000, pi1=0.4, pi0=0.02, k=0.5,
#'                 weight_se_f1=1, weight_se_prec=0, weight_se_rec=0)
#' ## example 3: now we want SE for F1 to be up to 3pp and SE for recall to be up to 5pp
#' ## and precision's SE is below 5pp. plus, assume in case multiple allocations verify this,
#' ## we want the allocation minimizing SE(F1)+0.5 SE(recall)+0.5 SE(precision).
#' ## assumptions about pi1, pi0 and k are as in example 2.
#' test_samplesize(se_f1=0.03, se_precision=0.05, se_recall=0.05,
#'                 max_N=4000, pi1=0.4, pi0=0.02, k=0.5,
#'                 weight_se_f1=1, weight_se_prec=0.5, weight_se_rec=0.5)
#' ## example 4: same as example 2 but we do not have a guess for pi0. however, we
#' ## know the classifier had recall of 0.6 in a dataset where 20% of observations
#' ## were positive
#' test_samplesize(se_f1=0.03, max_N=4000, pi1=0.4, k=0.5,
#'                 recall=0.6, external_positive_share=0.2,
#'                 weight_se_f1=1, weight_se_prec=0, weight_se_rec=0)
#' ## example 5: same as example 4 but we do not know the imbalance on the dataset
#' ## where recall was estimated, so we assume it to be same as in our evaluation
#' ## dataset
#' test_samplesize(se_f1=0.03, max_N=4000, pi1=0.4, k=0.5, recall=0.6,
#'                 weight_se_f1=1, weight_se_prec=0, weight_se_rec=0)
#' ## example 6: same as example 1 but we only care about precision and want
#' ## to get its SE below 3pp. then we do not need to specify pi0, recall,
#' ## external_k or external_positive_share
#' test_samplesize(se_precision=0.03, max_N=4000, pi1=0.4, k=0.5,
#'                 weight_se_f1=0, weight_se_prec=1, weight_se_rec=0)
#' @export
test_samplesize <- function(se_f1=NULL, se_precision=NULL, se_recall=NULL, max_N=NULL,
                            min_N=1, by_N=1, pi1=NULL, pi0=NULL, recall=NULL, k=NULL,
                            positive_share=NULL, external_k=NULL, external_positive_share=NULL,
                            weight_se_f1=1, weight_se_prec=1, weight_se_rec=1){
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

  if(!is.numeric(se_precision) & !is.numeric(se_recall) & !is.numeric(se_f1)){
    stop("No minimum standard error indicated")
  } else if(is.numeric(se_precision) & !is.numeric(se_recall) & !is.numeric(se_f1)){
    pi0 <- NA
    if(!is.numeric(weight_se_f1) | !is.numeric(pi1) | !is.numeric(k) |
       !is.numeric(max_N) | !is.numeric(min_N) | !is.numeric(by_N)){
      stop("Some parameters have not been rightly specified")
    }
    if(pi1<=0 | pi1>=1 | ((max_N %% 1) !=0) | ((min_N %% 1) !=0) | ((by_N %% 1) !=0)){
      stop("Some parameters have not been rightly specified")
    }
  } else {
    if(!is.numeric(weight_se_f1) | !is.numeric(weight_se_rec) | !is.numeric(weight_se_prec) | !is.numeric(pi1) |
       !is.numeric(pi0) | !is.numeric(k) | !is.numeric(max_N) | !is.numeric(min_N) | !is.numeric(by_N)){
      stop("Some parameters have not been rightly specified")
    }
    if(pi1<=0 | pi1>=1 | pi0<=0 | pi0>=1 | ((max_N %% 1) !=0) | ((min_N %% 1) !=0) | ((by_N %% 1) !=0)){
      stop("Some parameters have not been rightly specified")
    }
  }

  if(!is.numeric(se_f1)){
    weight_se_f1 <- 0
    se_f1 <- NA
  } else if(se_f1<=0 | se_f1>=1){
    stop("Invalid se_f1")
  }
  if(!is.numeric(se_recall)){
    weight_se_rec <- 0
    se_recall <- NA
  } else if(se_recall<=0 | se_recall>=1){
    stop("Invalid se_recall")
  }
  if(!is.numeric(se_precision)){
    weight_se_prec <- 0
    se_precision <- NA
  } else if(se_precision<=0|se_precision>=1){
    stop("Invalid se_precision")
  }

  recall <- k*pi1/(k*pi1+pi0)
  f1 <- 2*(recall*pi1)/(recall+pi1)
  f_star <- f1/(2-f1)

  precision_var_srs <- (pi1*(1-pi1)/((k/(k+1))*(1:max_N)))
  recall_var_srs <- (recall*(1-recall)/(((pi0+pi1*k)/(k+1))*(1:max_N)))
  f1_var_srs <- (f1*(1-f1)*(2-f1)^2)/(2*(1:max_N)*((pi0+k)/(k+1)))
  output_df_srs <- data.frame(precision_se_srs=sqrt(precision_var_srs),
                              recall_se_srs=sqrt(recall_var_srs),
                              f1_se_srs=sqrt(f1_var_srs),
                              sample_size_srs=1:max_N)
  # output_df_srs <- output_df_srs %>%
  #   dplyr::filter((is.na(se_precision) | (precision_se_srs < se_precision)) &
  #                 (is.na(se_recall) | (recall_se_srs < se_recall)) &
  #                 (is.na(se_f1) | (f1_se_srs < se_f1))) %>%
  #   dplyr::arrange(sample_size_srs) %>%
  #   dplyr::slice(1)
  output_df_srs <- output_df_srs[(is.na(se_precision) | output_df_srs$precision_se_srs < se_precision) &
                                 (is.na(se_recall) | output_df_srs$recall_se_srs < se_recall) &
                                 (is.na(se_f1) | output_df_srs$f1_se_srs < se_f1),]
  output_df_srs <- output_df_srs[order(output_df_srs$sample_size_srs), ]
  output_df_srs <- output_df_srs[1, , drop = FALSE]
  rownames(output_df_srs) <- 1

  for(n_sample in seq(min_N, max_N, by_N)){
    precision_var_strat <- (pi1*(1-pi1)/(1:n_sample))
    recall_var_strat <- (pi0/(k*pi1*((1+(pi0/(k*pi1)))^2)))^2*(((1-pi1)/(pi1*(1:n_sample))) + ((1-pi0)/(pi0*(n_sample-(1:n_sample)))))
    f1_var_strat <- (4/(1+f_star)^4)*f_star^2*(((1-pi1)/(pi1*(1:n_sample))) + ((1-pi0)*pi0/((pi0+k)^2*(n_sample-(1:n_sample)))))
    output_df_strat <- data.frame(precision_se_strat=sqrt(precision_var_strat),
                                  recall_se_strat=sqrt(recall_var_strat),
                                  f1_se_strat=sqrt(f1_var_strat),
                                  n_positives=1:n_sample,
                                  sample_size_strat=n_sample)
    # output_df_strat <- output_df_strat %>%
    #   dplyr::filter((is.na(se_precision) | (precision_se_strat < se_precision)) &
    #                   (is.na(se_recall) | (recall_se_strat < se_recall)) &
    #                   (is.na(se_f1) | (f1_se_strat < se_f1))) %>%
    #   dplyr::mutate(det = weight_se_prec*precision_se_strat + weight_se_rec*recall_se_strat + weight_se_f1*f1_se_strat) %>%
    #   dplyr::arrange(det) %>%
    #   dplyr::slice(1) %>%
    #   dplyr::select(-det)
    output_df_strat <- output_df_strat[(is.na(se_precision) | output_df_strat$precision_se_strat < se_precision) &
                                   (is.na(se_recall) | output_df_strat$recall_se_strat < se_recall) &
                                   (is.na(se_f1) | output_df_strat$f1_se_strat < se_f1),]
    output_df_strat$det <- weight_se_prec*output_df_strat$precision_se_strat + weight_se_rec*output_df_strat$recall_se_strat + weight_se_f1*output_df_strat$f1_se_strat
    output_df_strat <- output_df_strat[order(output_df_strat$det), ]
    output_df_strat <- output_df_strat[1, , drop = FALSE]
    output_df_strat$det <- NULL
    rownames(output_df_strat) <- 1
    if(nrow(output_df_strat)>0){
      break
    }
  }

  if(nrow(output_df_srs)==0){
    output_df_srs <- "Under SRS, cannot achieve desired SEs with so few observations: need to increase maximum sample size"
  }
  if(nrow(output_df_strat)==0){
    output_df_strat <- "Under stratified sampling, cannot achieve desired SEs with so few observations: need to increase maximum sample size"
  }

  output <- list(srs=output_df_srs, stratified=output_df_strat)
  return(output)
}


get_k <- function(positive_share){
  return(positive_share/(1-positive_share))
}

get_pi0 <- function(precision, recall, imbalance){
  return(((1/recall)-1)*imbalance*precision)
}

