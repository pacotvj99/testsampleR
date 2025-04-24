#' Estimate optimal allocation
#'
#' This function shows how many observations to sample per stratum under efficient stratified sampling
#'
#' This function can be used to determine how many observations to sample per stratum
#' under efficient stratified sampling. To do so, you should stratify by the
#' predicted probability of exhibiting the outcome, as produced by the classifier to be
#' evaluated: this should be the `stratifying` variable. The strata in `strata` should be
#' based on this variable, and neither the `stratifying` variable nor the `strata`
#' variable should contain NAs. At minimum, the `strata` variable should have one
#' stratum for the positive and one for the negatives. You can create additional
#' strata by binning the predicted probability, but importantly each stratum should
#' contain only positive or only negative observations. The function first determines
#' how many positive observations should be sampled in the test set, given parameter
#' values and objective function, calling \code{\link{optimal_n_positives}} internally.
#' This is the value that minimizes an objective function subject to a sample size
#' constraint (see documentation of function \code{\link{optimal_n_positives}} for details).
#' Then, the function performs proportional allocation within the positive and negative
#' subsamples separately, adjusting the number of sampled observations per bin so
#' that the overall sample size is as close to the target as possible.
#' The following parameters need to be entered to receive an output: (i) pi1 (precision,
#' TP/(TP+FP)), (ii) pi0 (FN/(FN+TN)) or recall (FN/(FN+TN)). These parameters should
#' be guessed or estimated on other data (pi0, pi1, recall). Note that if recall is entered instead
#' of pi0, then it is advised to also enter the imbalance ratio or positive share
#' of the dataset on which the recall guess is based (external_k or external_positive_share):
#' if this is not done, the imbalance in the test data is estimated by the function and used instead.
#'
#' @param data A data.frame from which to sample observations. One column should
#' contain the strata used for sampling.
#' @param N_sample A numeric value capturing the desired sample size.
#' @param strata A character value capturing the column name of the sampling strata.
#' Column should contain factor or character values. The strata should be based on
#' the same variable entered in `stratifying`.
#' @param stratifying A character value capturing the column name of the stratifying
#' variable. This column should contain continuous numeric values between 0 and 1,
#' and its values should reflect the predicted probability of exhibiting the outcome.
#' @param min_per_bin A numeric value capturing the minimum number of observations per bin.
#' Default is 1.
#' @param threshold A numeric value capturing the threshold for binary classification.
#' @param pi1 A numeric value capturing the expected precision (share of cases
#' exhibiting the outcome out of predicted positives). This value must be entered
#' for the function to produce an output.
#' @param pi0 A numeric value capturing the inverse of the expected precision for
#' the negative class (share of cases exhibiting the outcome out of predicted
#' negatives). If left blank, the function computes `pi0` from recall and imbalance.
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
#' @param weight_se_f1 A positive numeric value capturing how much to weigh
#' the SE for the F1 score in the objective function to be minimized. The default
#' for `weight_se_f1` is `1`.
#' @param weight_se_rec A positive numeric value capturing how much to weigh
#' the SE for recall in the objective function to be minimized. The default for
#' `weight_se_rec` is `0`.
#' @param weight_se_prec A positive numeric value capturing how much
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
#' @references
#' Tomas-Valiente, F. (2025). Uncertain performance: How to quantify uncertainty and
#' draw test sets when evaluating classifiers.
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
    stop("NAs in stratifying variable")
  }
  if(any(is.na(data$strata))){
    stop("NAs in strata variable")
  }
  data$stratifying <- unlist(data[,stratifying])
  data$strata <- unlist(data[,strata])

  if(threshold>1 | threshold<0){
    stop("Invalid threshold")
  }
  if(!is.numeric(N_sample) || !is.numeric(min_per_bin) || is.na(N_sample) || is.na(min_per_bin)){
    stop("Invalid binning specifications")
  } else if(((N_sample%%1)!=0) || ((min_per_bin%%1)!=0)){
    stop("Invalid binning specifications")
  }
  if(N_sample>nrow(data)){
    stop("Error: sample size is bigger than population data")
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
    stop("Some stratum contains both positive and negative observations")
  }

  #estimate optimal number of positives
  k <- sum(data$stratifying>=threshold)/sum(data$stratifying<threshold)
  optimal_n <- optimal_n_positives(N_sample, pi1, pi0, recall, k=k, positive_share=NULL,
                                   external_k, external_positive_share,
                                   weight_se_f1, weight_se_rec, weight_se_prec)
  if((!is.numeric(optimal_n))|(length(optimal_n)>1)){
    stop("Invalid optimal number of positives to be sampled")
  }

  #if optimally sample fewer positives than there are bins, coerce the sampled number of positives to the number of bins same for negatives
  if(optimal_n < length(unique(sam_pos$strata[!is.na(sam_pos$strata)]))){
    optimal_n <- length(unique(sam_pos$strata))*min_per_bin
    message("More positive bins than positives to be optimally sampled")
  } else if((N_sample - optimal_n) < length(unique(sam_neg$strata[!is.na(sam_neg$strata)]))){
    optimal_n <- N_sample - length(unique(sam_neg$strata))*min_per_bin
    message("More negative bins than positives to be optimally sampled")
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
#' The following parameters need to be entered to receive an output: (i) pi1 (precision,
#' TP/(TP+FP)), (ii) pi0 (FN/(FN+TN)) or recall (FN/(FN+TN)), (iii) imbalance ratio k
#' ((TP+FP)/(TN+FN)) or positive share ((TP+FP)/(TN+FN+TP+FP)). These parameters should
#' be guessed or estimated on other data (pi0, pi1, recall), or based on the predicted
#' labels in the population (k, positive share). Note that if recall is entered instead
#' of pi0, then it is advised to also enter the imbalance ratio or positive share
#' of the dataset on which the recall guess is based (external_k or external_positive_share):
#' if this is not done, the imbalance in the test data (k or positive_share) is used instead.
#'
#' @param N_sample A numeric value capturing the total test set sample size considered
#' as binding constraint (includes both positives and negatives).
#' @param pi1 A numeric value capturing the expected precision (share of cases
#' exhibiting the outcome out of predicted positives). This value must be entered
#' for the function to produce an output.
#' @param pi0 A numeric value capturing the inverse of the expected precision for
#' the negative class (share of cases exhibiting the outcome out of predicted
#' negatives). If left blank, the function computes `pi0` from recall and imbalance.
#' @param recall A numeric value capturing the expected value of recall (share of
#' positive cases among those exhibiting the outcome). This value only needs to be
#' entered if `pi0` is left blank.
#' @param k A numeric value capturing the imbalance ratio (number of predicted positives
#' over predicted negatives) in the population from which the test set should be
#' drawn. Argument `k` can be left blank, but then `positive_share` must be specified.
#' Typically, it should be below 1, as positives are defined as the rare class.
#' @param positive_share A numeric value capturing the positive share (number of
#' predicted positives out of all observations) in the population from which the
#' test set should be drawn. Argument `positive_share` can be left blank, but
#' then `k` must specified. Typically, it should be below 0.5, as positives are defined as the rare class.
#' @param external_k A numeric value capturing the imbalance ratio (number of predicted
#' positives over predicted negatives) in the population on which the recall estimate
#' is based. This only needs to be entered if `recall` is based on a population
#' that has a different imbalance than the test's population. The argument is
#' only required if `pi0` is not entered directly. If left blank, it will be
#' calculated internally using `external_positive_share`, if this is entered, or
#' assumed to be equal to the test population's imbalance.
#' Typically, it should be below 1, as positives are defined as the rare class.
#' @param external_positive_share A numeric value capturing the positive share
#' (number of predicted positives out of all observations) in the population on
#' which the recall estimate is based. This only needs to be entered if `recall`
#' is based on a population that has a different imbalance than the test's population.
#' The argument is only required if `pi0` is not entered directly. If left blank,
#' it will be calculated internally using `external_positive_share`, if this is
#' entered, or assumed to be equal to the test population's imbalance.
#' Typically, it should be below 0.5, as positives are defined as the rare class.
#' @param weight_se_f1 A positive numeric value capturing how much to weigh
#' the SE for the F1 score in the objective function to be minimized by the calculator.
#' The default for `weight_se_f1` is `1`.
#' @param weight_se_rec A positive numeric value capturing how much to weigh
#' the SE for recall in the objective function to be minimized by the calculator.
#' The default for `weight_se_rec` is `0`.
#' @param weight_se_prec A positive numeric value capturing how much
#' to weigh the SE for precision in the objective function to be minimized by the
#' calculator. The default for `weight_se_prec` is `0`.
#' @return A numeric value capturing how many positive observations should be sampled
#' given the parameter values and the objective function specified by the user.
#' @examples
#' ## example 1: how many out of 500 sampled obs should be positive if our goal is
#' ## minimizing the SE of F1? The answer here assumes that we expect precision of
#' ## 0.4, pi0 of 0.05, and 20% of positives
#' optimal_n_positives(N_sample = 500, pi1=0.4, pi0=0.05, positive_share=0.2,
#'                    weight_se_f1=1, weight_se_rec=0, weight_se_prec=0)
#' ## example 2: same but now we assume recall of 0.6, and that this was estimated on a
#' ## sample with 0.4 of positives
#' optimal_n_positives(N_sample = 500, pi1=0.4, recall=0.6, positive_share=0.2,
#'                    external_positive_share=0.2, weight_se_f1=1, weight_se_rec=0, weight_se_prec=0)
#' ## example 3: same as in example 1, but in the sample there is 1 positive for
#' ## every 2 negatives (k=1/2)
#' optimal_n_positives(N_sample = 500, pi1=0.4, recall=0.6, k=0.5, external_positive_share=0.2,
#'                     weight_se_f1=1, weight_se_rec=0, weight_se_prec=0)
#' ## example 4: same as in example 1, but we aim to minimize SE(F1)+0.5 SE(recall)+0.5 precision
#' optimal_n_positives(N_sample = 500, pi1=0.4, pi0=0.05, positive_share=0.2,
#'                     external_positive_share=0.2, weight_se_f1=1,
#'                     weight_se_rec=0.5, weight_se_prec=0.5)
#'
#' @export
optimal_n_positives <- function(N_sample=NULL, pi1=NULL, pi0=NULL, recall=NULL, k=NULL,
                                positive_share=NULL, external_k=NULL, external_positive_share=NULL,
                                weight_se_f1=1, weight_se_rec=0, weight_se_prec=0){
  #turn NaN and other weird into NA
  if(is.null(N_sample) || is.na(N_sample)){
    N_sample <- NA
  }
  if(is.null(pi1) || is.na(pi1)){
    pi1 <- NA
  }
  if(is.null(pi0) || is.na(pi0)){
    pi0 <- NA
  }
  if(is.null(recall) || is.na(recall)){
    recall <- NA
  }
  if(is.null(k) || is.na(k)){
    k <- NA
  }
  if(is.null(positive_share) || is.na(positive_share)){
    positive_share <- NA
  }
  if(is.null(external_k) || is.na(external_k)){
    external_k <- NA
  }
  if(is.null(external_positive_share) || is.na(external_positive_share)){
    external_positive_share <- NA
  }
  if(is.null(weight_se_f1) || is.na(weight_se_f1)){
    weight_se_f1 <- NA
  }
  if(is.null(weight_se_rec) || is.na(weight_se_rec)){
    weight_se_rec <- NA
  }
  if(is.null(weight_se_prec) || is.na(weight_se_prec)){
    weight_se_prec <- NA
  }

  #check that the imbalance is rightly specified
  if(!is.numeric(k) & is.numeric(positive_share)){
    if(positive_share<=0 | positive_share>=1){
      stop("Invalid positive_share: positive_share is above 1 or below 0")
    } else {
      if(positive_share>0.5){
        warning("positive_share above 0.5: you said there are more positives than negatives, but typically positives are defined as the rare class")
      }
      k <- get_k(positive_share)
    }
  } else if(is.numeric(k) & is.numeric(positive_share)){
    message("Both k and positive share specified: k used for analysis")
    if(k>1){
      warning("k above 1: you said there are more positives than negatives, but typically positives are defined as the rare class")
    }
  } else if(!is.numeric(k) & !is.numeric(positive_share)){
    stop("No information on imbalance: both k and positive_share are missing or invalid")
  } else if(k>1){
    warning("k above 1: you said there are more positives than negatives, but typically positives are defined as the rare class")
  }

  #checks on weights and sample sizes
  if(!is.numeric(weight_se_f1) | !is.numeric(weight_se_rec) | !is.numeric(weight_se_prec)){
    stop("Some weights have not been rightly specified: check weight_se_f1, weight_se_rec, weight_se_prec")
  } else if(weight_se_f1<0 | weight_se_rec<0 | weight_se_prec<0){
    stop("Some weights are negative: check weight_se_f1, weight_se_rec, weight_se_prec")
  }
  if(!is.numeric(N_sample)){
    stop("Missing sample size: check N_sample")
  } else if(((N_sample %% 1) !=0)){
    stop("Non-integer sample size: check N_sample")
  }

  #check pi1
  if(!is.numeric(pi1)){
    stop("Missing precision: pi1 is missing or invalid")
  } else if(pi1<=0 | pi1>=1){
    stop("Invalid precision: pi1 is above 1 or below 0")
  }

  #checks on pi0 and reconstruct if not entered or non-numeric
  if(is.numeric(pi0) & !is.na(pi0)){
    if(is.numeric(recall) & !is.na(recall)){
      message("Both pi0 and recall specified: only pi0 used for analysis")
    }
    if(pi0<=0 | pi0>=1){
      stop("Invalid pi0: pi0 is above 1 or below 0")
    }
  } else {
    if(!is.numeric(recall) | is.na(recall)){
      stop("Recall and pi0 were both missing or invalid: enter recall or pi0")
    } else {
      if(recall<=0 | recall >=1){
        stop("Invalid recall: recall is above 1 or below 0")
      } else {
        # first get the external_k
        if(!is.numeric(external_k) & is.numeric(external_positive_share)){
          if(external_positive_share<=0 | external_positive_share>=1){
            stop("Invalid external_positive_share: external_positive_share is above 1 or below 0")
          } else {
            if(external_positive_share>0.5){
              warning("external_positive_share above 0.5: you said there are more positives than negatives, but typically positives are defined as the rare class")
            }
            external_k <- get_k(external_positive_share)
          }
        } else if(is.numeric(external_k) & is.numeric(external_positive_share)){
          message("Both external_k and external_positive_share specified: external_k used for analysis")
        } else if(!is.numeric(external_k) & !is.numeric(external_positive_share)){
          message("In-sample imbalance assumed to estimate pi0")
          external_k <- k
        }
        if(external_k>1){
          warning("external_k above 1: you said there are more positives than negatives, but typically positives are defined as the rare class")
        }

        #then estimate pi0
        pi0 <- get_pi0(pi1, recall, external_k)
        if(!is.numeric(pi0) || is.na(pi0)){
          message("pi0 was calculated but is invalid")
        } else if(pi0<=0 | pi0>=1){
          stop("pi0 was calculated but is invalid")
        }
      }
    }
  }

  #hard check of key parameters
  if(!is.numeric(pi1) | !is.numeric(pi0) | !is.numeric(k)){
    stop("Some parameters have not been rightly specified")
  } else if(is.na(pi1) | is.na(pi0) | is.na(k)){
    stop("Some parameters have not been rightly specified")
  } else if(pi1<=0 | pi1>=1 | pi0<=0 | pi0>=1){
    stop("Some parameters have not been rightly specified")
  }

  #estimate optimal number of positives
  recall <- k*pi1/(k*pi1+pi0)
  if(recall<0.02){#check if very low recall implicitly assumed
    message("Very small recall implied for test set: consider revising assumptions")
  }
  f1 <- 2*(recall*pi1)/(recall+pi1)
  f_star <- f1/(2-f1)

  f1_variance <- (4/(1+f_star)^4)*f_star^2*(((1-pi1)/(pi1*(1:(N_sample-1)))) + ((1-pi0)*pi0/((pi0+k)^2*(N_sample-(1:(N_sample-1))))))
  recall_variance <- (pi0/(k*pi1*((1+(pi0/(k*pi1)))^2)))^2*(((1-pi1)/(pi1*(1:(N_sample-1)))) + ((1-pi0)/(pi0*(N_sample-(1:(N_sample-1))))))
  precision_variance <- (pi1*(1-pi1)/((1:(N_sample-1))))
  optimal_n <- which.min(weight_se_f1*sqrt(f1_variance) + weight_se_rec*sqrt(recall_variance) + weight_se_prec*sqrt(precision_variance))

  return(optimal_n)
}
