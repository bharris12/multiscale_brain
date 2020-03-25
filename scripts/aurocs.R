
#' predictors is a matrix where each column is a predictor and each row is a sample.
#' label_matrix is a binary matrix where columns are labels and each row is a sample.
#' 1 indicates that the sample on this row belongs to the label on this column.
compute_aurocs = function(predictors, label_matrix, return_tie_correction = FALSE) {
  predictors = as.matrix(predictors)
  label_matrix = as.matrix(label_matrix)
  n_positives = colSums(label_matrix)
  n_negatives = nrow(label_matrix) - n_positives
  ranks = matrixStats::colRanks(predictors, ties.method = "average", preserveShape=TRUE)
  sum_of_positive_ranks = crossprod(label_matrix, ranks)
  if (return_tie_correction) {
    tie_correction = compute_tie_correction(ranks)
  }
  colnames(sum_of_positive_ranks) = colnames(predictors)
  result = (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
  if (return_tie_correction) {
    return(list(aurocs = result, tie_corrections = tie_correction))
  } else {
    return(result)
  }
}

# For the following two functions, see
#
#   https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Normal_approximation_and_tie_correction
#
# The tie correction effectively computes lost variance because of ties (compared to discrete uniform).
# Computing the Wikipedia formula naively is slow, this method is equivalent and fast.
compute_tie_correction = function(ranks) {
  ranks = as.matrix(ranks)
  observed_var = matrixStats::colVars(ranks)
  max_var = var(seq_len(nrow(ranks)))
  return((max_var-observed_var) * 12 / nrow(ranks))
}

auroc_p_value = function(aurocs, label_matrix, two_tailed = TRUE, tie_correction = 0, log.p = FALSE) {
  p = colSums(label_matrix)
  n = nrow(label_matrix) - p
  
  # Careful: NAs arise from tie_correction (predictor with 0 variance)
  if (length(tie_correction) > 1) {
    Z = (aurocs - 0.5) * sqrt(12*n*p)
    Z = t(t(Z) / sqrt(n+p+1-tie_correction))
  } else {
    Z = (aurocs - 0.5) / sqrt((n+p+1-tie_correction)/(12*n*p))
  }
  
  result = Z
  if (two_tailed) {
    result[Z<=0] = pnorm(Z[Z<=0], log.p = log.p) * 2
    result[Z>0] = pnorm(Z[Z>0], lower.tail=FALSE, log.p = log.p) * 2
  } else {
    result = pnorm(Z, lower.tail=FALSE, log.p = log.p)
  }
  return(result)
}

design_matrix = function(cell_type, scale_columns=FALSE) {
  factors = levels(as.factor(cell_type))
  if (length(factors) > 1) {
    result = model.matrix(~as.factor(cell_type)-1)
  } else {
    result = matrix(1, nrow = length(cell_type), ncol = 1)
  }
  colnames(result) = factors
  if (scale_columns) {
    result = scale(result, center = FALSE, scale = colSums(result))
  }
  return(result)
}
