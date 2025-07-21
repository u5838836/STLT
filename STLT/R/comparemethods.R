#' @title Compare prediction performance of different models
#'
#' @description Compare prediction accuracy of qx between different methods using root mean squared error, mean absolute error, weighted root mean squared error and weighted mean absolute error.
#'
#' @param qx Vector of observed mortality rates; ages should be integer values and sequential.
#'
#' @param wx Vector of weights for the weighted prediction error measures. Usually the cohort size at each x. Can be a scalar indicating the cohort size at the first age - subsequent weights will then be calculated using the observed qx.
#'
#' @param pred_qx A matrix of predicted qx, corresponding to the ages of the observed qx. Each column should correspond to a different method used to predict qx.
#'
#' @return A matrix of prediction accuracy measures by model.
#'
#' @examples
#' qx <- c(0.1, 0.2, 0.15, 0.3)
#' pred_qx <- matrix(c(0.12, 0.11,0.21, 0.19,0.13, 0.18,0.31, 0.29), ncol = 2)
#' colnames(pred_qx) <- c("Method_A", "Method_B")
#' compare_methods(qx, pred_qx, wx = 10000)
#'
#' @export compare_methods

compare_methods<-function(qx, pred_qx, wx=NULL){
  if (length(qx) != nrow(pred_qx)) stop("qx length must match number of rows in pred_qx")

  # If no weights provided at all, use 1
  if (is.null(wx)) {
    wx <- 1
  }

  # If wx is a single number, generate weights using survival rule
  if (!is.null(wx) && length(wx) == 1) {
    w_full <- numeric(length(qx))
    w_full[1] <- wx
    for (i in 2:length(qx)) {
      w_full[i] <- w_full[i - 1] * (1 - qx[i - 1])
    }
    wx <- w_full
  }

  # Normalize weights to sum to 1
  w_norm <- wx / sum(wx)

  # Initialize results
  methods <- colnames(pred_qx)
  if (is.null(methods)) methods <- paste0("Method_", seq_len(ncol(pred_qx)))

  results <- data.frame(
    Method = methods,
    MAE = NA_real_,
    RMSE = NA_real_,
    Weighted_MAE = NA_real_,
    Weighted_RMSE = NA_real_
  )

  # Loop over prediction methods
  for (i in seq_len(ncol(pred_qx))) {
    pred <- pred_qx[, i]
    errors <- pred - qx
    abs_errors <- abs(errors)
    sq_errors <- errors^2

    results$MAE[i] <- mean(abs_errors)
    results$RMSE[i] <- sqrt(mean(sq_errors))
    results$Weighted_MAE[i] <- sum(w_norm * abs_errors)
    results$Weighted_RMSE[i] <- sqrt(sum(w_norm * sq_errors))
  }

  return(results)
}


