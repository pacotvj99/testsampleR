#' Implement stratified sampling
#'
#' This function implements stratified sampling under constant, proportional or optimal allocation
#'
#' This function can be used to implement the sampling process. It take as input the
#' dataset from which to sample and produces a smaller sampled dataset as output, i.e.
#' the test set. The function supports stratified sampling with constant allocation
#' across bins, proportional allocation, optimal two-bin allocation (with the
#' optimal number of positives for minimizing some linear combination of the SEs of
#' precision, recall and F1 score), and custom allocation to be manually inputed.
#' The function also allows users to flexibly select binning strategy, by specifying
#' how many bins they want to each side of the predicted probability threshold.
#' If users want to stratify on the basis of some variable other than the predicted
#' probability, they can construct the bins manually before using the function and
#' then enter a character variable as `stratifying` (this supports all types of
#' allocation except optimal). The output dataset contains the same columns as the
#' input plus: sampling indices, strata bins, and sampling probabilities which will
#' be used after annotation for constructing stratified estimates of performance metrics.
#' For optimal allocation, the following parameters need to be entered: (i) pi1 (precision,
#' TP/(TP+FP)), (ii) pi0 (FN/(FN+TN)) or recall (FN/(FN+TN)). These parameters should
#' be guessed or estimated on other data (pi0, pi1, recall). Note that if recall is entered instead
#' of pi0, then it is advised to also enter the imbalance ratio or positive share
#' of the dataset on which the recall guess is based (external_k or external_positive_share):
#' if this is not done, the imbalance in the test data is estimated by the function and used instead.
#'
#' @param data A data.frame from which to sample observations. One column should
#' contain the stratifying variable used to construct the strata for sampling.
#' @param stratifying A character value capturing the column name of the stratifying
#' variable. This column in `data` should contain continuous numeric values between
#' 0 and 1 (predicted probability of exhibiting the outcome), or else character labels
#' corresponding to the strata. If the variable captures predicted probabilities,
#' the positive class should typically be defined as the rare class. If user wants to
#' rely on custom strata, enter the name of character variable containing the strata.
#' If `allocation` is 'optimal', the variable indicated in `stratifying` should contain
#' continuous outcome probability: it should not be categorical.
#' @param N_sample A numeric value capturing the desired sample size. Not needed if
#' allocation is set to 'manual'.
#' @param na.rm A character value indicating how to deal with NA values in the
#' stratifying variable (indicated by `stratifying`). Options include: 'drop'
#' (this drops observations with NA in the stratifying variable), 'stop' (this
#' ensures the function stops running if there are any NAs in the stratifying variable),
#' 'impute' (this median-imputes the NAs in the stratifying variable) and encode
#' (this turns the NAs into a stratum of their own; this option is not supported
#' when using optimal). Default is 'stop'. When 'manual' allocation is selected,
#' users should create a stratum for NAs manually before calling the function.
#' @param allocation A character value capturing the type of allocation regime used
#' for stratified sampling. Options include: 'constant', 'proportional', 'optimal'
#' and 'manual'. If set to 'constant', the function calls \code{\link{constant_allocation}},
#' and implements constant allocation, aiming for all strata to get the same number
#' of observations in the test set. If set to 'proportional', the function calls
#' \code{\link{proportional_allocation}}, and implements proportional allocation,
#' aiming for each stratum to get a proportional number of observations in the test
#' set relative to the stratum's number of observations in the population. If set to
#' 'optimal', the function attempts to call \code{\link{optimal_allocation}} (only
#' if parameter values support this), and implements efficient two-bin allocation,
#' determining how many positives and negatives to optimally sample (or taking this
#' from the user-provided `n_positive`) and then performing proportional allocation
#' in each stratum. If set to 'manual', the function samples as many observations per
#' bin as indicated in the `manual_allocation` argument.
#' @param threshold A numeric value capturing the threshold for binary classification.
#' Not needed if allocation is set to 'manual'.
#' @param N_bins_left A numeric value capturing how many bins there should between
#' zero and `threshold` (closed and open respectively). Not needed if allocation is
#' set to 'manual'.
#' @param N_bins_right A numeric value capturing how many bins there should between
#' `threshold` and 1 (both closed). Not needed if allocation is set to 'manual'.
#' @param seed A numeric value that sets seed internally for sampling.
#' @param min_per_bin A numeric value capturing the minimum number of observations per bin.
#' Default is 1. Not needed if allocation is set to 'manual'.
#' @param manual_allocation A named vector, whose elements are numeric, and whose
#' names correspond to unique strata. This vector indicates how many observations to
#' sample, allowing the user to specify a specific custom stratified allocation.
#' This argument should only be provided when argument `allocation` is set to 'manual'.
#' @param pi1 A numeric value capturing the expected precision (share of cases
#' exhibiting the outcome out of predicted positives). This value is then used to
#' determine what the optimal allocation is. This argument only must be entered if
#' `allocation` is set to 'optimal' and no `caret` model object is provided in
#' `trained_model`.
#' @param pi0 A numeric value capturing the inverse of the expected precision for
#' the negative class (share of cases exhibiting the outcome out of predicted
#' negatives). This value is then used to determine what the optimal allocation is.
#' This argument only must be entered if `allocation` is set to 'optimal' and no
#' `caret` model object is provided in `trained_model` and `recall` is not entered.
#' @param trained_model A `caret` model object, used to estimate `pi1` and `pi0` if
#' these are not manually provided, and if simultaneously `allocation` is set to
#' 'manual'. The model can only be used if cross validation was used during
#' training: then, excluded folds are used to guess `pi1` and `pi0`.
#' @param recall A numeric value capturing the expected value of recall (share of
#' positive cases among those exhibiting the outcome). This argument only needs to
#' be entered if `allocation` is set to 'optimal' and `pi0` is not entered and no
#' `caret` model object is provided in `trained_model`
#' @param external_k A numeric value capturing the imbalance ratio (number of predicted
#' positives over predicted negatives) in the population on which the recall estimate
#' is based. This only needs to be entered if `recall` is based on a population
#' that has a different imbalance than the test's population. The argument is
#' only required if `pi0` is not entered directly and no `caret` model object is
#' provided in `trained_model` and `external_positive_share` is not entered, while
#' allocation is set to 'optimal'. If defined, it should typically be below 1, as
#' positives are defined as the rare class.
#' @param external_positive_share A numeric value capturing the positive share
#' (number of predicted positives out of all observations) in the population on
#' which the recall estimate is based. This only needs to be entered if `recall`
#' is based on a population that has a different imbalance than the test's population.
#' The argument is only required if `pi0` is not entered directly and no `caret`
#' model object is provided in `trained_model` and `external_k` is not entered,
#' while allocation is set to 'optimal'. If defined, it should typically be below 0.5, as
#' positives are defined as the rare class.
#' @param n_positive A numeric value indicating how many positive observations to
#' include in the stratified sample. Positives and negatives are then sampled
#' separately using stratified sampling with proportional allocation in the positive
#' and negative subsamples respectively. Note `n_positive` should be an integer,
#' and it should only be entered if `allocation` is set to 'optimal'.
#' @param weight_se_f1 A positive numeric value capturing how much to weigh
#' the SE for the F1 score in the objective function to be minimized for optimal
#' allocation (see details in \code{\link{optimal_allocation}}). The default
#' for `weight_se_f1` is `1`. This should only be entered if `allocation` is set
#' to 'optimal'.
#' @param weight_se_prec A positive numeric value capturing how much to
#' weigh the SE for precision in the objective function to be minimized for
#' optimal allocation (see details in \code{\link{optimal_allocation}}). The
#' default for `weight_se_prec` is `0`. This should only be entered if `allocation`
#' is set to 'optimal'.
#' @param weight_se_rec A positive numeric value capturing how much to weigh
#' the SE for recall in the objective function to be minimized for optimal
#' allocation (see details in \code{\link{optimal_allocation}}). The default for
#' `weight_se_rec` is `0`. This should only be entered if `allocation` is set to
#' 'optimal'.
#' @return A data.frame object that contains the test set sample to be annotated.
#' This data.frame contains same columns as the original data.frame inputted in the
#' `data` argument, but contains fewer rows since it is just a sample. Additionally it
#' contains the following column: `stratifying` (numeric variable with the values of the
#' stratifying variable if the input corresponded to a numeric stratifying variable;
#' character variable indicating the strata if this was the input), `id_sampling`
#' (variable indicating which rows were included in the sample), `strata` (character
#' variable indicating what stratum each sampled observation belongs to), and `Prob`
#' (numeric variable indicating the sampling probability for any given observation in
#' the test sample).
#' @examples
#' ## example 1: draw a test set of 1000 observations, by stratified sampling on
#' ## predicted probability bins, with 5 bins left of 0.5 threshold and 5 bins right
#' ## of this threshold with constant allocation. The function also supports different
#' ## thresholds and different number of bins to either side of cutoff (including
#' ## asymmetric bin numbers), but 5-5 is a reasonable default for small samples:
#' ## if bigger sample, can have more bins. threshold should be the value at which
#' ## one starts to deem an observation positive (here 0.5, but alternative values
#' ## are possible).
#' data(pop_df)
#' out <- testsampler(data=pop_df, stratifying='score', N_sample=1000, allocation='constant',
#'                    threshold=0.5, N_bins_left=5, N_bins_right=5)
#' ## example 2: draw a test set of 1000 observations, by stratified sampling on
#' ## predicted probability bins, with 5 bins left of 0.5 threshold and 5 bins right
#' ## of this threshold with proportional allocation.
#' data(pop_df)
#' out <- testsampler(data=pop_df, stratifying='score', N_sample=1000, allocation='proportional',
#'                    threshold=0.5, N_bins_left=5, N_bins_right=5)
#' ## example 3: draw a test set of 1000 observations, by stratified sampling on
#' ## predicted probability bins, with 5 bins left of 0.5 threshold and 5 bins right
#' ## of this threshold with optimal allocation. first, optimal number of positives and negatives
#' ## is determined, and then proportional allocation within positives/negatives. assume we're
#' ## minimizing the objective function SE(F1) + 0.5 SE(precision) + 0.5 SE(recall).
#' ## as assumptions, we have guesses for pi1 (TP/TP+FP) and pi0 (FN/FN+TN) of 0.4 and 0.02.
#' data(pop_df)
#' out <- testsampler(data=pop_df, stratifying='score', N_sample=1000, threshold=0.5,
#'                    N_bins_left=5, N_bins_right=5, pi1=0.4, pi0=0.2, weight_se_f1=1,
#'                    weight_se_rec=0.5, weight_se_prec=0.5)
#' ## example 4: like example 3, but instead we only care about the SE of the f1.
#' ## Plus, we want simple two-bin positive-negative stratification (0.5 threshold).
#' data(pop_df)
#' out <- testsampler(data=pop_df, stratifying='score', N_sample=1000, threshold=0.5,
#'                    N_bins_left=1, N_bins_right=1, pi1=0.4, pi0=0.2, weight_se_f1=1,
#'                    weight_se_rec=0, weight_se_prec=0)
#' ## example 5: like example 3, but assuming we know recall was 0.6 in some other
#' ## dataset where the imbalance ratio had 1 positive for every 2 negatives so
#' ## external_k is 1/2=0.5
#' data(pop_df)
#' out <- testsampler(data=pop_df, stratifying='score', N_sample=1000, threshold=0.5,
#'                    N_bins_left=5, N_bins_right=5, pi1=0.4, recall=0.6, external_k=0.5,
#'                    weight_se_f1=1, weight_se_rec=0.5, weight_se_prec=0.5)
#' ## example 6: we decide manually how many positives and negatives to sample (500
#' ## each) and then sample with proportional allocation with the positive and negative
#' ## strata (5 bins in each).
#' data(pop_df)
#' out <- testsampler(data=pop_df, stratifying='score', N_sample=1000, threshold=0.5,
#'                    N_bins_left=5, N_bins_right=5, n_positive=500)
#' ## example 7: like example 1, we use constant allocation stratification, but on
#' ## some non-numeric categorical variable. therefore, no threshold needed
#' data(pop_df)
#' pop_df$somevar <- as.character(rbinom(nrow(pop_df), 2, 0.4)) #create random stratifying variable
#' out <- testsampler(data=pop_df, stratifying='somevar', allocation='constant',
#'                    N_sample=1000, N_bins_left=5, N_bins_right=5)
#' ## example 8: like example 7, we stratify the sampling on some custom categorical
#' ## variable, but specifying a custom allocation
#' data(pop_df)
#' pop_df$somevar <- as.character(rbinom(nrow(pop_df), 2, 0.4)) #create stratifying variable
#' custom <- c(200, 200, 600)
#' names(custom) <- c("0", "1", "2") # these names need to match the values of somevar
#' out <- testsampler(data=pop_df, stratifying='somevar', allocation='manual',
#'                    manual_allocation=custom)
#' ## example 9: optimal stratification with categorical labels is only possible if
#' ## each stratum contains only positives or only negatives. This might happen when
#' ## you want to create custom bins, e.g. bins that are not to a constant distance
#' ## to each other. You can implement this using manual allocation in the following
#' ## way, assuming pi1=0.4 and pi0=0.2.
#' data(pop_df)
#' # stratify in custom bins
#' pop_df$custom_bins <- cut(pop_df$score, breaks=c(0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1))
#' # determine optimal allocation
#' custom <- optimal_allocation(data = pop_df, N_sample = 500, strata = 'custom_bins',
#'                               stratifying = 'score', min_per_bin=1, threshold=0.5,
#'                               pi1=0.4, pi0=0.2, weight_se_f1=1, weight_se_rec=0,
#'                               weight_se_prec=0)
#' # draw test set with manual allocation
#' out <- testsampler(data=pop_df, stratifying='custom_bins', allocation='manual',
#'                    manual_allocation=custom)
#' @references
#' Tomas-Valiente, F. (2025). Uncertain performance: How to quantify uncertainty and
#' draw test sets when evaluating classifiers.
#' @export
testsampler <- function(data, stratifying, N_sample=NA, na.rm="stop", allocation="optimal",
                        threshold=0.5, N_bins_left=5, N_bins_right=5, seed=1234, min_per_bin=1,
                        manual_allocation=NULL, pi1=NULL, pi0=NULL, trained_model=NULL,
                        recall=NULL, external_k=NULL, external_positive_share=NULL,
                        n_positive=NULL, weight_se_f1=1, weight_se_rec=0, weight_se_prec=0){
  if("data.frame" %in% class(data)){
    class(data) <- "data.frame"
  } else {
    stop("Invalid input data")
  }
  if(!is.character(stratifying)){
    stop("Invalid stratifying argument")
  } else if(sum(colnames(data)==stratifying)!=1){
    stop("Column with stratifying variable cannot be uniquely identified")
  }

  #ensure variable is called stratifying
  for(conflicting_column in c("Prob", "strata", "stratifying", "id_sampling")){
    if(conflicting_column!="stratifying" | stratifying!="stratifying"){
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
  data$stratifying <- unlist(data[,stratifying])

  #check other arguments
  if(!is.null(n_positive) & allocation!="optimal"){
    message("Number of positives specified, so allocation should be set to 'optimal'")
    allocation <- "optimal"
  }
  if(!(allocation %in% c("optimal", "proportional", "constant", "manual"))){
    stop("Invalid allocation argument")
  } else if(allocation!="manual"){
    if(!is.numeric(threshold) || threshold>1 || threshold<0 || is.na(threshold)){
      stop("Invalid threshold")
    }
    if(!is.numeric(N_sample) || is.na(N_sample) || ((N_sample%%1)!=0)){
      stop("Invalid specifications: N_sample missing or invalid")
    }
    if(!is.numeric(N_bins_left) || is.na(N_bins_left) || ((N_bins_left%%1)!=0)){
      stop("Invalid specifications: N_bins_left missing or invalid")
    }
    if(!is.numeric(N_bins_right) || is.na(N_bins_right) || ((N_bins_right%%1)!=0)){
      stop("Invalid specifications: N_bins_right missing or invalid")
    }
    if(!is.numeric(min_per_bin) || is.na(min_per_bin) || ((min_per_bin%%1)!=0)){
      stop("Invalid specifications: min_per_bin missing or invalid")
    }
    if(is.numeric(manual_allocation)){
      stop("Specific manual allocation was entered, but allocation type was not manual")
    }
  }

  if(!is.null(n_positive)){
    n_positive <- as.numeric(n_positive)
    if((is.na(n_positive)) | (n_positive %% 1 != 0)){
      stop("Number of positives specified, but not integer")
    }
    if(n_positive>N_sample){
      stop("Number of positives specified is greater than sample size")
    }
    if(!is.null(trained_model) | is.numeric(pi1) | is.numeric(pi0) | is.numeric(recall)){
      message("Additional parameters and number of positives are both specified: number of positives used for allocation")
    }
  }

  #drop NAs or impute them
  if(!(na.rm %in% c("drop", "stop", "impute", "encode"))){
    stop("Invalid na.rm argument")
  }
  if(any(is.na(data$stratifying))){
    if(allocation=="manual"){
      stop("NAs in stratifying variable are not supported under manual allocation: create a stratum for them or use another allocation")
    } else if(na.rm=="drop") {
      message("NAs detected and removed from the data.")
      data <- data[!is.na(data$stratifying),]
    } else if(na.rm=="impute") {
      if(is.numeric(data$stratifying)){
        message("NAs detected and imputed using median")
        data$stratifying <- ifelse(is.na(data$stratifying), stats::median(data$stratifying, na.rm=TRUE), data$stratifying)
      } else {
        stop("NAs in stratifying variable that cannot be median-imputed!")
      }
    } else if(na.rm=="stop") {
      stop("NAs in stratifying variable!")
    } else if(na.rm=="encode" & allocation == "optimal") {
      stop("NAs in stratifying variable cannot be encoded under optimal stratification")
    }
  }

  #create indices for sampling
  data$id_sampling <- as.character(1:nrow(data))

  # perform stratification
  if(na.rm=="encode"){
    missing_data <- data[is.na(data$stratifying),]
    data <- data[!is.na(data$stratifying),]
  }
  if(is.numeric(data$stratifying)){
    if(min(data$stratifying)<0 | max(data$stratifying)>1){
      stop("Stratifying variable is numeric but not probability")
    }

    # stratify by positives and negatives separately
    if(any(data$stratifying<threshold) & any(data$stratifying>=threshold)){
      sam_neg <- data[data$stratifying<threshold,]
      sam_neg$strata <- cut(sam_neg$stratifying, breaks = seq(0,threshold,threshold/N_bins_left), include.lowest = TRUE)
      sam_pos <- data[data$stratifying>=threshold,]
      sam_pos$strata <- cut(sam_pos$stratifying, breaks = seq(threshold,1,(1-threshold)/N_bins_right), include.lowest = TRUE)
      data <- rbind(sam_neg, sam_pos)
      if(nrow(sam_neg)<nrow(sam_pos)){
        warning("There are more positives than negatives in the data: typically, positives are defined as the rare class")
      }
    } else if(!any(data$stratifying<threshold) & any(data$stratifying>=threshold)){
      warning("There are only positives in the data: check the threshold is right")
      data$strata <- cut(data$stratifying, breaks = seq(threshold,1,(1-threshold)/N_bins_right), include.lowest = TRUE)
    } else if(any(data$stratifying<threshold) & !any(data$stratifying>=threshold)){
      warning("There are only negatives in the data: check the threshold is right")
      data$strata <- cut(data$stratifying, breaks = seq(0,threshold,threshold/N_bins_left), include.lowest = TRUE)
    }
    # clean strata content
    data$strata <- gsub("\\[", "\\(", gsub("\\]", "\\)", data$strata))
  } else {
    data$strata <- as.character(data$stratifying)
    if(allocation=="optimal"){
      stop("Non-numeric stratifying variable cannot sustain optimal allocation.")
    }
  }
  if(na.rm=="encode"){
    if(nrow(missing_data)>0){
      message("NAs in stratifying variable have been encoded as a stratum")
      missing_data$strata <- "Recoded_missings"
      data <- rbind(data, missing_data)
    }
  }
  data <- data[order(data$strata), ]

  # perform the allocation
  if(allocation == "constant"){
    myallocation <- constant_allocation(data=data, N_sample=N_sample, strata="strata", min_per_bin=min_per_bin)
  } else if(allocation == "proportional"){
    myallocation <- proportional_allocation(data=data, N_sample=N_sample, strata="strata", min_per_bin=min_per_bin)
  } else if(allocation == "optimal"){
    if(!is.numeric(n_positive)){
      pi1 <- ifelse(is.null(pi1) || is.na(pi1), NA, pi1)
      pi0 <- ifelse(is.null(pi0) || is.na(pi1), NA, pi0)
      if(is.numeric(pi0)){
        if(is.numeric(recall)){
          message("Both pi0 and recall specified: only pi0 used for analysis")
        }
      }
      #if input a model rather than pi1 and pi0 estimates, then estimate pi1 and pi0 inside the function
      if((is.na(pi1) | is.na(pi0)) & ("train" %in% class(trained_model))){
        try({
          pi1 <- as.numeric(mean(trained_model$resampledCM$cell4/(trained_model$resampledCM$cell2+trained_model$resampledCM$cell4), na.rm=T))
          pi0 <- as.numeric(mean(trained_model$resampledCM$cell3/(trained_model$resampledCM$cell3+trained_model$resampledCM$cell1), na.rm=T))
        }, silent = TRUE)
        pi1 <- ifelse(is.null(pi1) || is.na(pi1), NA, pi1)
        pi0 <- ifelse(is.null(pi0) || is.na(pi1), NA, pi0)
      }
    }

    if(is.numeric(n_positive)){
      message("Optimal allocation with pre-specified number of positives")
      myallocation_negatives <- proportional_allocation(data=data[data$stratifying<threshold,], N_sample=(N_sample-n_positive), strata="strata", min_per_bin=min_per_bin)
      myallocation_positives <- proportional_allocation(data=data[data$stratifying>=threshold,], N_sample=n_positive, strata="strata", min_per_bin=min_per_bin)
      myallocation <- c(myallocation_negatives, myallocation_positives)
    } else if(max(data$stratifying, na.rm=T)<threshold | min(data$stratifying, na.rm=T)>=threshold){
      warning("Optimal allocation could not be performed because all observations fall to the same side of threshold: resorting to proportional allocation")
      myallocation <- proportional_allocation(data=data, N_sample=N_sample, strata="strata", min_per_bin=min_per_bin)
    } else if(pi1>=1 | pi1<=0 | is.na(pi1)){
      stop("Optimal allocation could not be performed because pi1 was not entered or is invalid: reenter pi1 or choose another allocation")
    } else if(!is.numeric(pi0) || is.na(pi0)){
      if(is.numeric(recall) && !is.na(recall)){
        if(recall>=1 | recall<=0){
          stop("Optimal allocation could not be performed because recall was invalid: reenter recall or choose another allocation")
        } else {
          myallocation <- optimal_allocation(data = data, N_sample = N_sample, strata = "strata", stratifying = "stratifying",
                                             min_per_bin = min_per_bin, threshold = threshold,
                                             pi1 = pi1, pi0 = pi0, recall = recall, external_k = external_k,
                                             external_positive_share = external_positive_share,
                                             weight_se_f1 = weight_se_f1, weight_se_rec = weight_se_rec,
                                             weight_se_prec = weight_se_prec)
        }
      } else {
        stop("Optimal allocation could not be performed because neither pi0 nor recall was not entered: enter pi0/recall or choose another allocation")
      }
    } else if(pi0>=1 | pi0<=0){
      stop("Optimal allocation could not be performed because pi0 was invalid: reenter pi0 or choose another allocation")
    } else {
      myallocation <- optimal_allocation(data=data, N_sample=N_sample, strata="strata", stratifying="stratifying",
                                         min_per_bin=min_per_bin, threshold=threshold,
                                         pi1=pi1, pi0=pi0, recall=recall, external_k=external_k,
                                         external_positive_share=external_positive_share, weight_se_f1=weight_se_f1,
                                         weight_se_rec=weight_se_rec, weight_se_prec=weight_se_prec)
    }
  } else if(allocation == "manual"){
    if(is.numeric(N_sample) & !is.na(N_sample)){
      warning("Allocation was manual, but sample size was entered: sample size ignored")
    }
    myallocation <- manual_allocation
    if(length(unique(data$strata)) != length(myallocation)){
      stop("Manual allocation does not contain as many entries as strata")
    }
    if(length(unique(data$strata)) != length(names(myallocation))){
      stop("Manual allocation does not contain as many names as strata")
    }
    if(length(unique(names(myallocation))) != length(names(myallocation))){
      stop("Manual allocation contains non-unique names")
    }
    if(!is.numeric(myallocation)){
      stop("Non-numeric manual allocation ")
    }
    if(any(is.na(myallocation))){
      stop("Some missings in manual allocation ")
    }
    myallocation <- myallocation[order(names(myallocation))]
  }

  #actual sampling
  strata_count <- table(data$strata)
  strata_count <- strata_count[order(names(strata_count))]
  myallocation <- myallocation[order(names(myallocation))]
  if(any(is.na(myallocation))){
    stop("Missings in allocation")
  }
  if(any(names(strata_count)!=names(myallocation))){
    stop("The strata in the data do not match the strata in the manual allocation")
  }
  if(length(unique(names(myallocation)))!=length(strata_count)){
    stop("Some stratum was entered multiple times")
  }
  if((any(myallocation<min_per_bin & strata_count>=min_per_bin)  & allocation!="manual") ||
     any((myallocation<1) | (strata_count < myallocation))){
    stop("Invalid allocation")
  }

  which_in_sample <- c()
  #sample as many as allocated for each stratum and generate sampling weights
  data$Prob <- NA
  for(stratum_i in names(myallocation)){
    set.seed(seed)
    which_in_sample <- c(which_in_sample, sample(data$id_sampling[data$strata==stratum_i], myallocation[stratum_i]))
    data$Prob[data$strata==stratum_i] <- myallocation[stratum_i]/strata_count[stratum_i]
  }
  if(any(is.na(data$Prob))){
    stop("Failure of common support")
  }

  design_str <- data[data$id_sampling %in% which_in_sample,]
  set.seed(seed)
  design_str <- design_str[sample(1:nrow(design_str), nrow(design_str), replace = FALSE),]
  return(design_str)
}
