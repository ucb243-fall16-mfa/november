#' @title eigentbl
#' @description Return table with eigens, cumulative, % inertia, and cumulative % of inertia
#' @param mfa an R object
#' @export
eigentbl <- function(mfa) {
  sing_vals <- round(sqrt(mfa$eigenvalues),3)
  eigen <- round(mfa$eigenvalues,3)
  eigen_cum <- round(cumsum(mfa$eigenvalues),3)
  perc_inertia <- round(mfa$eigenvalues / sum(mfa$eigenvalues),2) * 100
  perc_cum <- round(cumsum(mfa$eigenvalues / sum(mfa$eigenvalues)),2) * 100
  df <- cbind(sing_vals,eigen,eigen_cum,perc_inertia,perc_cum)
  return(df)
}

#' @title ob2dim
#' @description Compute the contribution of an observation to a dimension
#' @param mfa an R object
#' @export
ob2dim <- function(mfa) {
  M <- sapply(1:length(mfa$eigenvalues), function(i) (1/12 * mfa$factor_scores[,i]^2) / mfa$eigenvalues[i])
  return(M)
}

#' @title var2dim
#' @description Compute the contribution of a variable to a dimension
#' @param mfa an R object
#' @export
var2dim <- function(mfa) {
  mfa$A %*% t(mfa$loadings^2)
}

#' @title tbl2dim
#' @description Compute the contribution of a table to a dimension
#' @param mfa an R object
#' @export
tbl2dim <- function(mfa) {
  # call to var_to_dim helper function
  v2d <- var2dim(mfa)
  # create a matrix to store the results
  M <- matrix(nrow=length(mfa$sets),ncol=12)
  for (i in 1:length(mfa$sets)) {
    M[i,] <- apply(v2d,2,function(x) sum(x[mfa$sets[[i]]]))
  }
  return(M)
}

#' @title RV
#' @description Compute the similarity between two tables using the Rv coefficient
#' @param tbl1 a matrix or data frame
#' @param tbl2 a matrix or data frame
#' @export
RV <- function(tbl1, tbl2) {
  # encode the tables as matrices
  tbl1 <- as.matrix(tbl1)
  tbl2 <- as.matrix(tbl2)
  num <- sum(diag((tcrossprod(tbl1) %*% tcrossprod(tbl2))))
  den <- sqrt(sum(diag(tcrossprod(tbl1) %*% tcrossprod(tbl1))) *
              sum(diag(tcrossprod(tbl2) %*% tcrossprod(tbl2))))
  return(num/den)
}

#' @title RVtbl
#' @description Compute the similarity matrix between multiple sets of tables
#' @param data data set (matrix or data frame)
#' @param sets list of vectors indicating the sets of variables
#' @export
RVtbl <- function(data, sets) {
  # encode the data as a matrix
  X <- sapply(sets, function(x) data[,x][])
  # merge list of matrices into a single matrix
  X <- do.call(cbind, X)
  # specify that data is formatted as a matrix
  X <- as.matrix(X)
  # create a matrix to store the results
  M <- diag(length(sets))
  # call helper function to adjust the indices of sets
  sets.len <- sapply(sets, function(x) length(x))
  sets <- adjust_sets(sets, sets.len)
  # use the "comb" base function to generate all possible combinations of tables
  combos <- combn(1:length(sets),2)
  for (i in 1:ncol(combos)) {
    a <- combos[1,i]
    b <- combos[2,i]
    M[a,b] <- RV(X[,sets[[a]]],X[,sets[[b]]])
    M[b,a] <- RV(X[,sets[[a]]],X[,sets[[b]]])
  }
  return(M)
}
