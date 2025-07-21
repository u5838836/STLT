#' @title Evaluate prediction performance
#'
#' @description Evaluate prediction performance of any time-varying model using root mean squared error, mean absolute error, weighted root mean squared error and weighted mean absolute error.
#'
#' @param qx Matrix of observed mortality rates
#'
#' @param wx Matrix or vector of weights for the weighted prediction error measures. If vector, vector should be the starting cohort size for each cohort.
#'
#' @param pred_qx List of matrices (one matrix for each method) of predicted qx with rows corresponding to the ages of the observed qx. Each column should correspond to a different cohort.
#'
#' @return A vector of prediction accuracy measures.
#'
#' @examples
#' qx <- matrix(runif(15, 0.01, 0.3), nrow = 5, ncol = 3)
#' pred1 <- qx + matrix(rnorm(15, 0, 0.02), nrow = 5)
#' pred2 <- qx + matrix(rnorm(15, 0, 0.03), nrow = 5)
#'
#' # Starting cohort sizes
#' start_weights <- c(1000, 1500, 2000)
#'
#' # Evaluate
#' error_measures(
#'   qx = qx,
#'   pred_qx = list(Method_A = pred1, Method_B = pred2),
#'   wx = start_weights
#' )
#'
#' @export error_measures

error_measures <- function(qx, pred_qx, wx = NULL) {
  # qx: matrix (rows = age, columns = cohorts)
  # pred_qx: list of matrices (same dim as qx), one per method
  # wx: either a matrix (same dim), or a vector of starting cohort sizes

  # Input checks
  if (!is.matrix(qx)) stop("qx must be a matrix")
  if (!is.list(pred_qx)) pred_qx <- list(Method_1 = pred_qx)

  n_age <- nrow(qx)
  n_cohort <- ncol(qx)
  n_methods <- length(pred_qx)

  # Generate weights
  if (is.null(wx)) {
    w_mat <- matrix(1, nrow = n_age, ncol = n_cohort)
  } else if (is.matrix(wx)) {
    if (!all(dim(wx) == dim(qx))) stop("Weight matrix must have same dimensions as qx")
    w_mat <- wx
  } else if (length(wx) == n_cohort) {
    # Derive weights per cohort
    w_mat <- matrix(NA, nrow = n_age, ncol = n_cohort)
    for (j in seq_len(n_cohort)) {
      w_mat[1, j] <- wx[j]
      for (i in 2:n_age) {
        w_mat[i, j] <- w_mat[i - 1, j] * (1 - qx[i - 1, j])
      }
    }
  } else {
    stop("wx must be a matrix or a vector with length equal to number of cohorts")
  }

  # Normalize weights across all cells
  w_norm <- w_mat / sum(w_mat, na.rm = TRUE)

  # Results data frame
  results <- data.frame(
    Method = names(pred_qx),
    MAE = NA_real_,
    RMSE = NA_real_,
    Weighted_MAE = NA_real_,
    Weighted_RMSE = NA_real_
  )

  # Loop over methods
  for (m in seq_len(n_methods)) {
    pred <- pred_qx[[m]]

    if (!all(dim(pred) == dim(qx))) stop(paste("Prediction matrix for method", m, "has incorrect dimensions"))

    error <- pred - qx
    abs_err <- abs(error)
    sq_err <- error^2

    results$MAE[m] <- mean(abs_err, na.rm = TRUE)
    results$RMSE[m] <- sqrt(mean(sq_err, na.rm = TRUE))

    results$Weighted_MAE[m] <- sum(w_norm * abs_err, na.rm = TRUE)
    results$Weighted_RMSE[m] <- sqrt(sum(w_norm * sq_err, na.rm = TRUE))
  }

  return(results)
}

