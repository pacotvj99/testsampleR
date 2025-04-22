#' Estimate performance metrics from confusion matrix
#'
#' This function estimates performance metrics given a test set confusion matrix
#'
#' This function can be used to estimate metrics from a confusion matrix. The function
#' can be applied regardless of how the test set was sampled, provided the estimated
#' share of observation in each cell of the confusion matrix is unbiased:
#' raw shares are unbiased if test set was SRS, stratified estimators are needed if
#' test set was drawn with stratified sampling. Note that the sum of the proportion of
#' observations across all four cells must be one. The function is called internally
#' within \code{\link{stratified_metrics}}. The following metrics are supported:
#' precision, recall, F1 score for the positive class, precision of negative class,
#' recall of negative class, F1 score for the negative class, accuracy, MCC (Matthews
#' correlation coefficient), Cohen's kappa, (unweighted) macro-averaged F1 score,
#' weighted F1 score, BM (bookmaker's informedness). The output also shows the
#' proportions of true negatives, false negatives, false positives and true positives.
#'
#' @param tr0_pred0 A numeric value capturing the share of true negatives.
#' @param tr1_pred0 A numeric value capturing the share of false negatives.
#' @param tr0_pred1 A numeric value capturing the share of false positives.
#' @param tr1_pred1 A numeric value capturing the share of true positives.
#' @return A data.frame of performance metrics. The following
#' metrics are supported: precision, recall, F1 of positive class, accuracy, precision
#' of negative class, recall of negative class, F1 of negative class, MCC (Matthews
#' correlation coefficient), Cohen's kappa, (unweighted) macro-averaged F1 score,
#' weighted F1 score, BM (bookmaker's informedness). The output also shows the
#' proportions of true negatives, false negatives, false positives and true positives.
#' @seealso \code{\link{stratified_metrics}}
#' @examples
#' ## enter the estimates for the test set confusion matrix cells
#' ## get estimates of key performance metriccs
#' do_metrics(tr0_pred0 = 0.5, tr1_pred0 = 0.1, tr0_pred1 = 0.1, tr1_pred1 = 0.3)
#' @export
do_metrics <- function(tr0_pred0, tr1_pred0, tr0_pred1, tr1_pred1){
  if(is.na(tr0_pred0) | is.na(tr1_pred0) | is.na(tr0_pred1) | is.na(tr1_pred1)){
    stop("Invalid confusion matrix.")
  }
  if(!is.numeric(tr0_pred0) | !is.numeric(tr1_pred0) | !is.numeric(tr0_pred1) | !is.numeric(tr1_pred1)){
    stop("Invalid confusion matrix.")
  }
  if(tr0_pred0<0 | tr0_pred0>1 | tr1_pred0<0 | tr1_pred0>1 | tr0_pred1<0 | tr0_pred1>1 | tr1_pred1<0 | tr1_pred1>1){
    stop("Invalid confusion matrix.")
  }
  if(round(tr0_pred0+tr1_pred0+tr0_pred1+tr1_pred1, 12)!=1){
    stop("Confusion matrix cells do not add up to one.")
  }

  precision <- tr1_pred1/(tr1_pred1+tr0_pred1)
  recall <- tr1_pred1/(tr1_pred1+tr1_pred0)
  f1 <- ifelse(precision==0 & recall==0, 0, 2/((1/precision)+(1/recall)) )
  accuracy <- tr0_pred0+tr1_pred1
  MCC <- (tr1_pred1*tr0_pred0 - tr1_pred0*tr0_pred1)/sqrt((tr1_pred1+tr0_pred1)*(tr1_pred1+tr1_pred0)*(tr0_pred0+tr0_pred1)*(tr0_pred0+tr1_pred0))
  kappa <- 2*(tr1_pred1*tr0_pred0 - tr0_pred1*tr1_pred0)/((tr1_pred1+tr0_pred1)*(tr0_pred1+tr0_pred0) + (tr1_pred1+tr1_pred0)*(tr0_pred0+tr1_pred0))
  precision_rev <- tr0_pred0/(tr0_pred0+tr1_pred0)
  recall_rev <- tr0_pred0/(tr0_pred0+tr0_pred1)
  f1_neg <- ifelse(precision_rev==0 & recall_rev==0, 0, (2/((1/precision_rev)+(1/recall_rev))))
  f1_unweighted <- (f1+ f1_neg)/2
  f1_weighted <- ((tr1_pred1+tr1_pred0)*f1)+((tr0_pred0+tr0_pred1)*f1_neg)
  bm <- recall + (tr0_pred0/(tr0_pred0+tr0_pred1)) - 1

  output <- data.frame(f1=f1, precision=precision, recall=recall,
                       negative_precision=precision_rev, negative_recall=recall_rev,
                       negative_f1=f1_neg, accuracy=accuracy, MCC=MCC, Kappa=kappa,
                       f1_unweighted=f1_unweighted, f1_weighted=f1_weighted, BM=bm,
                       TN=tr0_pred0, FN=tr1_pred0, FP=tr0_pred1, TP=tr1_pred1)
  return(output)
}

#' Estimate performance metrics from test set
#'
#' This function estimates performance metrics given an annotated test set
#'
#' This function can be used on a test set drawn from simple random sampling, or from
#' stratified sampling if sampling probabilities are provided. The function unbiasedly
#' estimates the proportion of observations in each cell of the test set confusion
#' matrix, using stratified estimators if needed, and then produces key metrics based
#' on them by internally calling \code{\link{do_metrics}}. The input data should
#' contain a column capturing whether each observation exhibits the outcome (1/TRUE)
#' or not (0/FALSE), and a column capturing the predicted label (either predicted
#' probability, or binary label in either numeric or logical format). The following
#' metrics are supported: precision, recall, F1 of positive class, accuracy, precision
#' of negative class, recall of negative class, F1 of negative class, MCC (Matthews
#' correlation coefficient), Cohen's kappa, (unweighted) macro-averaged F1 score,
#' weighted F1 score, BM (bookmaker's informedness). The output also shows the
#' proportions of true negatives, false negatives, false positives and true positives.
#'
#' @param data A data.frame containing the annotated test set.
#' @param truth A character value capturing the column name of the annotated labels.
#' Column should capture whether each observation exhibits the outcome (1/TRUE)
#' or not (0/FALSE).
#' @param pred A character value capturing the column name of the classifier's
#' predictions. Column should capture each observation's predicted label (either
#' numeric or logical), or the predicted probability of exhibiting the outcome (numeric).
#' Default is NULL, which assumes equal sampling probability across observations.
#' @param probs A character value capturing the column name of the sampling probabilities.
#' @param threshold A numeric value capturing the threshold for binary classification.
#' @param se A logical value whether SEs should be calculated. If all
#' observations have the same sampling probability, the function calls `se_srs` to
#' estimate SEs using the analytical formula. In this case, Wilson CI lower and
#' upper bounds for precision, recall and F1 score are also outputted. If observations
#' do not have the same sampling probability, then the function calls `se_strat` to
#' compute stratified SEs: this assumes that the stratification was done on the basis
#' of positives and negatives (if proportional stratification was used within the
#' positives and negatives separately, the SEs are conservative). If stratification
#' was done on some other variable, then `se` should be set to FALSE, and bootstrapping
#' should be used instead.
#' @param bs A numeric value indicating the number of bootstrapping iterations used
#' to compute confidence intervals. If `NULL`, then bootstrapping is not implemented.
#' For stratified sampling, bootstrapping is stratified by stratum too.
#' @param strata A character value capturing the column name of the sampling strata.
#' Column should contain factor or character values. This only needs to be entered if
#' `bs` is not `NULL` and the sampling was stratified (sampling probabilities are not
#' constant across observations).
#' @param seed A numeric value used as seed for the bootstrapping.
#' @param alpha A numeric value used capturing the type-I error rate for constructing
#' confidence intervals (applies to both bootstrapping and analytic CIs).
#' @return A data.frame of performance metrics. The following
#' metrics are supported: precision, recall, F1 of positive class, accuracy, precision
#' of negative class, recall of negative class, F1 of negative class, MCC (Matthews
#' correlation coefficient), Cohen's kappa, (unweighted) macro-averaged F1 score,
#' weighted F1 score, BM (bookmaker's informedness). The output also shows the
#' proportions of true negatives, false negatives, false positives and true positives.
#' If `bs` was specified, then the resulting data.frame also produces the bootstrapped
#' SEs and CIs (columns ending in `boot_se`,`boot_lwr`, `boot_upr`): these are based
#' on stratified bootstrapping if the test set was drawn through stratified sampling.
#' If `se` was set to `TRUE`, analytic SEs are included in the output data.frame too
#' for precision, recall and F1 score (columns `f1_se`, `recall_se`, `precision_se`).
#' These SEs are valid for SRS and two-bin efficient stratification by positives vs
#' negatives. If proportional stratification was used when samplingpositives and
#' negatives, the analytic SEs are conservative. If other stratification regime was
#' used, these SEs should not be trusted and bootstrapped SEs/CIs should be considered
#' instead.
#' @examples
#' ## example 1: compute metrics after stratified sampling using probability weights
#' data(sample_df)
#' stratified_metrics(data = sample_df, truth = "truth", pred = "score", probs="Prob",
#'                    threshold = 0.5, strata = "strata", se=FALSE)
#' ## example 2: same as example 1 but with bootstrapped SEs/CIs
#' data(sample_df)
#' stratified_metrics(data = sample_df, truth = "truth", pred = "score", probs="Prob",
#'                    threshold = 0.5, strata = "strata", se=FALSE, bs=200)
#' ## example 3: same as example 1 but assuming the test set is a SRS from the
#' ## population. since it's a SRS, analytic SEs and CIs are valid and computed
#' data(sample_df)
#' sample_df$Prob <- 1
#' stratified_metrics(data = sample_df, truth = "truth", pred = "score", probs="Prob",
#'                    threshold = 0.5, strata = "strata", se=TRUE)
#' ## alternatively, can achieve the same with
#' stratified_metrics(data = sample_df, truth = "truth", pred = "score", probs=NULL,
#'                    threshold = 0.5, strata = "strata", se=TRUE)
#' ## can additionally get bootstrapped SEs and CIs with
#' stratified_metrics(data = sample_df, truth = "truth", pred = "score", probs=NULL,
#'                    threshold = 0.5, strata = "strata", se=TRUE, bs=200)
#' ## example 4: same as example 1 but computing SEs analytically.
#' ## this assumes that the test set is based on positive-negative two-bin
#' ## stratification. this provides conservative inference if the test set is drawn
#' ## by optimal stratification (first picking the number of positives and negatives
#' ## as in two-bin stratification, then proportional stratification within positive
#' ## and negative subsamples). the SEs are invalid if the stratification regime does
#' ## not stratify on the predicted probability of the outcome.
#' stratified_metrics(data = sample_df, truth = "truth", pred = "score", probs="Prob",
#'                    threshold = 0.5, strata = "strata", se=TRUE)

#' @seealso \code{\link{do_metrics}}, \code{\link{se_srs}}, \code{\link{se_strat}}
#' @export
stratified_metrics <- function(data, truth, pred, probs=NULL, threshold, se=TRUE,
                               bs=NULL, strata, seed=1234, alpha=0.05){
  if("data.frame" %in% class(data)){
    class(data) <- "data.frame"
  } else {
    stop("Invalid input data")
  }

  if(is.null(probs)){
    warning("No probabilities in data: constant sampling probability assumed. This is a problem if stratified sampling was used for drawing test set, but it is fine if the test set is a SRS")
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
  if(any(is.na(data$Prob))){
    stop("Some sampling probabilities are missing")
  }

  if(!(se %in% c(TRUE, FALSE))){
    stop("Invalid se argument")
  }
  if(!is.null(bs)){
    if(!is.numeric(seed)){
      stop("Invalid seed")
    }
    if(!is.numeric(bs) || is.na(bs)){
      stop("Invalid number of bootstrapping iterations")
    }
    if((bs %% 1)!=0){
      stop("Invalid number of bootstrapping iterations")
    }
    if(length(unique(data$Prob))>1){
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
      } else if(any(is.na(data$strata))){
        stop("Some strata are missing")
      }
    }
  }
  if(!is.numeric(alpha)){
    stop("Invalid alpha")
  } else if(alpha<=0 | alpha>=1){
    stop("Invalid alpha")
  }

  if(!is.character(truth)){
    stop("Invalid truth argument")
  }
  if(sum(colnames(data)==truth)!=1){
    stop("Column with annotated labels cannot be uniquely identified")
  }
  if(!is.character(pred)){
    stop("Invalid pred argument")
  }
  if(sum(colnames(data)==pred)!=1){
    stop("Column with predicted labels cannot be uniquely identified")
  }
  if(is.na(threshold) || !is.numeric(threshold) || threshold>1 || threshold<0){
    stop("Invalid threshold")
  }

  data$truth <- as.numeric(unlist(data[,truth]))
  data$var1 <- as.numeric(unlist(data[,pred]))
  if(any(!(data$truth %in% c(0,1)))){
    stop("Invalid annotated labels")
  }
  if(any(is.na(data$var1))){
    stop("Some predicted labels are missing")
  }
  if(any(data$var1<0 | data$var1>1)){
    stop("Invalid predicted labels")
  }

  # produce stratified estimates
  data$matrix_cell <- factor(ifelse(data$truth==0,
                                    ifelse(data$var1<threshold, "tr0_pred0", "tr0_pred1"),
                                    ifelse(data$var1<threshold, "tr1_pred0", "tr1_pred1")),
                             levels=c("tr0_pred0", "tr0_pred1", "tr1_pred0", "tr1_pred1"))
  # data <- data[!is.na(data$matrix_cell),]
  tr0_pred0 <- sum((data$matrix_cell == "tr0_pred0") * (1 / data$Prob)) / sum(1 / data$Prob)
  tr1_pred0 <- sum((data$matrix_cell == "tr1_pred0") * (1 / data$Prob)) / sum(1 / data$Prob)
  tr0_pred1 <- sum((data$matrix_cell == "tr0_pred1") * (1 / data$Prob)) / sum(1 / data$Prob)
  tr1_pred1 <- sum((data$matrix_cell == "tr1_pred1") * (1 / data$Prob)) / sum(1 / data$Prob)

  # analytic SEs
  output <- do_metrics(tr0_pred0, tr1_pred0, tr0_pred1, tr1_pred1)
  if(se){
    if(length(unique(data$Prob))==1){
      se_df <- se_srs(tr0_pred0, tr1_pred0, tr0_pred1, tr1_pred1, nrow(data))
      ci_df <- ci_srs_wilson(tr0_pred0, tr1_pred0, tr0_pred1, tr1_pred1, nrow(data), alpha=alpha)
      se_df <- cbind(se_df, ci_df)
    } else {
      message("Analytical SEs assume that stratification was done on the basis of positives vs negatives.\nFor other stratification regimes, use bootstrapping.")
      N_positive <- sum(data$var1>=threshold)
      N_negative <- sum(data$var1<threshold)
      se_df <- se_strat(tr0_pred0, tr1_pred0, tr0_pred1, tr1_pred1, N_positive, N_negative)
    }
    output <- cbind(output, se_df)
  }

  # bootstrapped SEs and CIs if request
  if(is.numeric(bs)){
    bs_output <- data.frame()
    set.seed(seed)
    for(j in 1:bs){
      if(length(unique(data$Prob))==1){
        which_sample <- sample(nrow(data), nrow(data), replace=TRUE)
        data_resample <- data[which_sample, ]
      } else {
        data_resample <- data.frame()
        for(stratum in unique(data$strata)){
          data_stratum <- data[data$strata==stratum,]
          which_sample <- sample(nrow(data_stratum), nrow(data_stratum), replace=TRUE)
          data_resample <- rbind(data_resample, data_stratum[which_sample, ])
        }
      }
      tr0_pred0 <- sum((data_resample$matrix_cell == "tr0_pred0") * (1 / data_resample$Prob)) / sum(1 / data_resample$Prob)
      tr1_pred0 <- sum((data_resample$matrix_cell == "tr1_pred0") * (1 / data_resample$Prob)) / sum(1 / data_resample$Prob)
      tr0_pred1 <- sum((data_resample$matrix_cell == "tr0_pred1") * (1 / data_resample$Prob)) / sum(1 / data_resample$Prob)
      tr1_pred1 <- sum((data_resample$matrix_cell == "tr1_pred1") * (1 / data_resample$Prob)) / sum(1 / data_resample$Prob)
      bs_output <- rbind(bs_output, do_metrics(tr0_pred0, tr1_pred0, tr0_pred1, tr1_pred1))
    }

    boot_se <- as.data.frame(t(apply(bs_output, 2, stats::sd)))
    colnames(boot_se) <- paste0(colnames(boot_se), "_boot_se")
    boot_lwr <- as.data.frame(t(apply(bs_output, 2, function(x) stats::quantile(x, probs = alpha / 2))))
    colnames(boot_lwr) <- paste0(colnames(boot_lwr), "_boot_lwr")
    boot_upr <- as.data.frame(t(apply(bs_output, 2, function(x) stats::quantile(x, probs = 1 - (alpha / 2)))))
    colnames(boot_upr) <- paste0(colnames(boot_upr), "_boot_upr")
    output <- cbind(output, boot_se, boot_lwr, boot_upr)
  }

  # reorder column of the dataset
  mycolnames <- colnames(output)[which(colnames(output)=="f1"):which(colnames(output)=="TP")]
  if(se){
    mycolnames <- c(sort(c(mycolnames[1:3], colnames(se_df))), mycolnames[4:length(mycolnames)])
  }
  if(is.numeric(bs)){
    bootnames <- colnames(output)[grepl("boot", colnames(output))]
    mycolnames2 <- unlist(lapply(mycolnames, function(somename){c(somename, bootnames[grepl(paste0("^", somename, "_boot"), bootnames)])}))
  } else {
    mycolnames2 <- mycolnames
  }
  if(!setequal(mycolnames2, colnames(output)) || length(unique(mycolnames2))!=ncol(output)){
    stop("Error in reordering final output")
  } else {
    output <- output[,mycolnames2]
  }

  return(output)
}
