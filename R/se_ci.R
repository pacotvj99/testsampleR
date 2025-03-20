############# Perhaps I should allow this to be done based on expected pi1, pi0, k like in test_samplesize
# se_srs_precision <- function(precision, positive_share, N_sample){
#   return(sqrt(precision*(1-precision)/(positive_share*N_sample)))
# }
# se_srs_recall <- function(recall, yes_share, N_sample){
#   return(sqrt(recall*(1-recall)/(yes_share*N_sample)))
# }
# se_srs_f1 <- function(f1, N_true){
#   return(sqrt(recall*(1-recall)/N_true))
# }

#' Estimate standard errors for simple random sample
#'
#' This function analytically estimates SEs from confusion matrix drawn by simple random sampling
#'
#' This function takes a confusion matrix and sample size as input, and produces
#' SEs with Wilson method for the F1 score, precision and recall. These SEs are
#' based on assuming the test set was drawn by simple random sampling: the CIs
#' should not be trusted if sampling was conducted in a different way.
#'
#' @param tr0_pred0 A numeric value capturing the share of true negatives.
#' @param tr1_pred0 A numeric value capturing the share of false negatives.
#' @param tr0_pred1 A numeric value capturing the share of false positives.
#' @param tr1_pred1 A numeric value capturing the share of true positives.
#' @param N_sample A numeric value capturing the total number size of the test set.
#' @return A data.frame containing the analytic SEs for the F1 score, precision and recall.
#' These CIs are based on assuming the test set was drawn by simple random sampling:
#' the CIs should not be trusted if sampling was conducted in a different way.
#' @examples
#' ## enter the estimates for the confusion matrix cells and sample size
#' ## get the SE analytically under assumption that test set is a SRS
#' se_srs(tr0_pred0 = 0.5, tr1_pred0 = 0.1, tr0_pred1 = 0.1, tr1_pred1 = 0.3, N_sample = 1000)
#' @export
se_srs <- function(tr0_pred0, tr1_pred0, tr0_pred1, tr1_pred1, N_sample){
  if(!is.numeric(tr0_pred0) | !is.numeric(tr1_pred0) | !is.numeric(tr0_pred1) | !is.numeric(tr1_pred1)){
    stop("Invalid confusion matrix.")
  }
  if(tr0_pred0<0 | tr0_pred0>1 | tr1_pred0<0 | tr1_pred0>1 | tr0_pred1<0 | tr0_pred1>1 | tr1_pred1<0 | tr1_pred1>1){
    stop("Invalid confusion matrix.")
  }
  if(round(tr0_pred0+tr1_pred0+tr0_pred1+tr1_pred1, 12)!=1){
    stop("Confusion matrix cells do not add up to one.")
  }
  if((N_sample %% 1)!=0){
    stop("Invalid sample size")
  }
  precision <- tr1_pred1/(tr1_pred1+tr0_pred1)
  recall <- tr1_pred1/(tr1_pred1+tr1_pred0)
  f1 <- ifelse(precision==0 & recall==0, 0, 2/((1/precision)+(1/recall)))

  precision_se <- sqrt(precision*(1-precision)/((tr0_pred1+tr1_pred1)*N_sample))
  recall_se <- sqrt(recall*(1-recall)/((tr1_pred0+tr1_pred1)*N_sample))
  f1_se <- sqrt(f1*(1-f1)*(2-f1)^2/(2*N_sample*(1-tr0_pred0)))

  output <- data.frame(f1_se=f1_se, precision_se=precision_se, recall_se=recall_se)
  return(output)
}

#' Estimate standard errors for stratified random sample
#'
#' This function analytically estimates SEs from confusion matrix drawn by simple random sampling
#'
#' This function takes a confusion matrix and sample size as input, and produces
#' SEs with Wilson method for the F1 score, precision and recall. These SEs are
#' based on assuming two-bin stratified random sampling, where the two bins are
#' positive and negative observations: the SEs should not be trusted if sampling
#' was conducted in a different way (either not by two-bin positives vs negatives
#' stratification, or by stratifying on some other variable).
#'
#' @param tr0_pred0 A numeric value capturing the share of true negatives.
#' @param tr1_pred0 A numeric value capturing the share of false negatives.
#' @param tr0_pred1 A numeric value capturing the share of false positives.
#' @param tr1_pred1 A numeric value capturing the share of true positives.
#' @param N_positive A numeric value capturing the number size of positive observations
#' sampled into the test set.
#' @param N_negative A numeric value capturing the number size of negative observations
#' sampled into the test set.
#' @return A data.frame containing the analytic SEs for the F1 score, precision and
#' recall under two-bin stratification of positives and negatives. These SEs are not
#' valid if the test set was not drawn through stratified sampling (e.g SRS) or if
#' it was drawn through some other stratification regime (either not by two-bin
#' positives vs negatives stratification, or by stratifying on some other variable).
#' @examples
#' ## enter the estimates for the test set confusion matrix cells and sample sizes
#' ## get the SE analytically under assumption that test set is a two-bin
#' ## stratified sample (positives and negatives)
#' se_strat(tr0_pred0 = 0.5, tr1_pred0 = 0.1, tr0_pred1 = 0.1, tr1_pred1 = 0.3,
#'          N_positive = 600, N_negative = 400)
#' @export
se_strat <- function(tr0_pred0, tr1_pred0, tr0_pred1, tr1_pred1, N_positive, N_negative){
  if(!is.numeric(tr0_pred0) | !is.numeric(tr1_pred0) | !is.numeric(tr0_pred1) | !is.numeric(tr1_pred1)){
    stop("Invalid confusion matrix.")
  }
  if(tr0_pred0<0 | tr0_pred0>1 | tr1_pred0<0 | tr1_pred0>1 | tr0_pred1<0 | tr0_pred1>1 | tr1_pred1<0 | tr1_pred1>1){
    stop("Invalid confusion matrix.")
  }
  if(round(tr0_pred0+tr1_pred0+tr0_pred1+tr1_pred1, 12)!=1){
    stop("Confusion matrix cells do not add up to one.")
  }
  if(((N_positive %% 1)!=0) || ((N_negative %% 1)!=0)){
    stop("Invalid sample size")
  }

  k <- (tr0_pred1 + tr1_pred1)/(tr0_pred0 + tr1_pred0)
  pi1 <- tr1_pred1/(tr1_pred1+tr0_pred1)
  pi0 <- tr1_pred0/(tr0_pred0+tr1_pred0)
  recall <- tr1_pred1/(tr1_pred1+tr1_pred0)
  f1 <- ifelse(pi1==0 & recall==0, 0, 2/((1/pi1)+(1/recall)))
  f_star <- f1/(2-f1)

  precision_se <- sqrt(pi1*(1-pi1)/N_positive)
  recall_se <- sqrt((pi0/(k*pi1*((1+(pi0/(k*pi1)))^2)))^2*(((1-pi1)/(pi1*N_positive)) + ((1-pi0)/(pi0*N_negative))))
  f1_se <- sqrt((4/(1+f_star)^4)*f_star^2*(((1-pi1)/(pi1*N_positive)) + ((1-pi0)*pi0/((pi0+k)^2*N_negative))))

  output <- data.frame(f1_se=f1_se, precision_se=precision_se, recall_se=recall_se)
  return(output)
}

#' Estimate Wilson confidence intervals for simple random sample
#'
#' This function estimates CIs with Wilson method from confusion matrix drawn by simple random sampling
#'
#' This function takes a confusion matrix and sample size as input, and produces
#' CIs with Wilson's score method for the F1 score, precision and recall. These CIs
#' are based on assuming the test set was drawn by simple random sampling: the CIs
#' should not be trusted if sampling was conducted in a different way.
#'
#' @param tr0_pred0 A numeric value capturing the share of true negatives.
#' @param tr1_pred0 A numeric value capturing the share of false negatives.
#' @param tr0_pred1 A numeric value capturing the share of false positives.
#' @param tr1_pred1 A numeric value capturing the share of true positives.
#' @param N_sample A numeric value capturing the total number size of the test set.
#' @param alpha A numeric value used capturing the type-I error rate for constructing
#' confidence intervals (applies to both bootstrapping and analytic CIs).
#' @return A data.frame containing the upper and lower endpoints of the Wilson CIs for
#' the F1 score, precision and recall. These CIs are based on assuming the test set
#' was drawn by simple random sampling: the CIs should not be trusted if sampling
#' was conducted in a different way.
#' @examples
#' ## enter the estimates for the test set confusion matrix cells and sample size
#' ## get the Wilson CIs analytically under assumption that test set is a SRS
#' ci_srs_wilson(tr0_pred0 = 0.5, tr1_pred0 = 0.1, tr0_pred1 = 0.1, tr1_pred1 = 0.3,
#' N_sample = 1000, alpha = 0.05)
#' @export
ci_srs_wilson <- function(tr0_pred0, tr1_pred0, tr0_pred1, tr1_pred1, N_sample, alpha=0.05){
  if(!is.numeric(tr0_pred0) | !is.numeric(tr1_pred0) | !is.numeric(tr0_pred1) | !is.numeric(tr1_pred1)){
    stop("Invalid confusion matrix.")
  }
  if(tr0_pred0<0 | tr0_pred0>1 | tr1_pred0<0 | tr1_pred0>1 | tr0_pred1<0 | tr0_pred1>1 | tr1_pred1<0 | tr1_pred1>1){
    stop("Invalid confusion matrix.")
  }
  if(round(tr0_pred0+tr1_pred0+tr0_pred1+tr1_pred1, 12)!=1){
    stop("Confusion matrix cells do not add up to one.")
  }
  if((N_sample %% 1)!=0){
    stop("Invalid sample size")
  }
  if(!is.numeric(alpha)){
    stop("Invalid alpha")
  } else if(alpha<=0 | alpha>=1){
    stop("Invalid alpha")
  }

  precision <- tr1_pred1/(tr1_pred1+tr0_pred1)
  recall <- tr1_pred1/(tr1_pred1+tr1_pred0)
  f1 <- ifelse(precision==0 & recall==0, 0, 2/((1/precision)+(1/recall)))
  f_star <- f1/(2-f1)
  n_pred_yes <- (tr0_pred1+tr1_pred1)*N_sample
  n_real_yes <- (tr1_pred0+tr1_pred1)*N_sample
  n_f1 <- (tr1_pred0+tr0_pred1+tr1_pred1)*N_sample

  precision_ci <- wilson_ci(precision, n_pred_yes, alpha=alpha)
  recall_ci <- wilson_ci(recall, n_real_yes, alpha=alpha)
  wilson_fstar_ci <- wilson_ci(f_star, n_f1, alpha=alpha)
  wilson_f1_upper <- 2*wilson_fstar_ci[2]/(1+wilson_fstar_ci[2])
  wilson_f1_lower <- 2*wilson_fstar_ci[1]/(1+wilson_fstar_ci[1])

  output <- data.frame(f1_wilson_lwr=wilson_f1_lower, f1_wilson_upr=wilson_f1_upper,
                       precision_wilson_lwr=precision_ci[1], precision_wilson_upr=precision_ci[2],
                       recall_wilson_lwr=recall_ci[1], recall_wilson_upr=recall_ci[2])
  return(output)
}

#' Estimate Wilson confidence intervals
#'
#' This helper function estimates CIs with Wilson method
#'
#' This function takes the estimate and sample size needed for the Wilson method
#' and computes the CI's endpoints by solving the associated second-order equation.
#'
#' @param estimate A numeric value capturing the point estimate.
#' @param sample_size_wilson A numeric value capturing the total number of observations
#' that should be considered for the Wilson CI. This should be the number of positives
#' if the estimate is precision, the number of observations that actually exhibit the
#' outcome if the estimate is recall, and the number of observations that are true positives
#' or false positives or false negatives if the estimate is F1.
#' @param alpha A numeric value used capturing the type-I error rate for constructing
#' confidence intervals (applies to both bootstrapping and analytic CIs).
#' @return A vector containing the lower and upper endpoints of the Wilson CI.
#' @examples
#' wilson_ci(estimate=0.4, sample_size_wilson=200, alpha=0.05)
#' @noRd
wilson_ci <- function(estimate, sample_size_wilson, alpha){
  if(!is.numeric(alpha)){
    stop("Invalid alpha")
  } else if(alpha<=0 | alpha>=1){
    stop("Invalid alpha")
  }
  if(!is.numeric(estimate) | !is.numeric(sample_size_wilson)){
    stop("Invalid input to Wilson CIs")
  }
  zstat <- stats::qnorm(1-(alpha/2))
  wilson_factor <- zstat^2/sample_size_wilson
  wilson_a <- (1+wilson_factor)
  wilson_b <- -(2*estimate+wilson_factor)
  wilson_c <- estimate^2
  wilson_upper <- ((-wilson_b)+sqrt(wilson_b^2-4*wilson_a*wilson_c))/(2*wilson_a)
  wilson_lower <- ((-wilson_b)-sqrt(wilson_b^2-4*wilson_a*wilson_c))/(2*wilson_a)
  return(c(wilson_lower, wilson_upper))
}

