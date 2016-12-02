#' @title mfa
#' @description Creates an object of class \code{"mfa"}
#' @param data data set (matrix or data frame)
#' @param sets list of vectors indicating the sets of variables
#' @param ncomps integer indicating how many components (i.e. factors) should be extracted
#' @param center either a logical value or a numeric vector of length equal to # of active vars
#' @param scale either a logical value or a numeric vector of length equal to # of active vars
#' @return an object of class coin
#' @export
#' @examples
#' # default
#' mfa1 <- mfa(data, sets=list(2:7,8:13,14:19), ncomps = NULL, center = TRUE, scale = TRUE))

mfa <- function(data, sets, ncomps, center = TRUE, scale = TRUE) {
  # call helper function to check that params are in proper format
  check_params(data, sets)

  # get descriptor vector
  names <- names(df)[unlist(sets)]
  # pull out data columns specified in sets
  X <- sapply(sets, function(x) data[,x][])
  # merge list of matrices into a single matrix
  X <- do.call(cbind, X)
  # specify that data is formatted as a matrix
  X <- as.matrix(X)

  # center and scale values
  if (center) X <- apply(X, 2, function(x) x - mean(x))
  if (scale) X <- apply(X, 2, function(x) x / sqrt(sum(x * x)))

  # function call to adjust sets to dimensions of reduced matrix
  sets.len <- sapply(sets, function(x) length(x))
  sets <- adjust_sets(sets, sets.len)

  # function call to compute the alpha weight vector
  a <- get_alphas(X,sets)
  # function call to convert alpha vector into A vector
  A <- get_A(a, sets.len)
  # compute M vector for diagonal mass matrix
  M <- diag(rep(1 / nrow(X), nrow(X)))

  # function to compute generalized PCA
  A.tilde <- sqrt(M) %*% X %*% sqrt(A)
  A.tilde.svd <- svd(A.tilde)
  U.tilde <- solve(sqrt(M)) %*% A.tilde.svd$u
  V.tilde <- solve(sqrt(A)) %*% A.tilde.svd$v
  D.tilde <- A.tilde.svd$d
  loadings <- t(V.tilde)
  factor_scores <- U.tilde %*% diag(D.tilde)

  # function call to compute partial factor score array
  P <- get_P(X, sets, a, V.tilde)

  # create the MFA object
  object <- list(
    eigenvalues = D.tilde^2,
    factor_scores = factor_scores,
    partial_factor_scores = P,
    loadings = loadings,
    A = A,
    labels = df[,1],
    sets = sets,
    names = names
  )

  # set the class to "mfa"
  class(object) <- "mfa"
  return(object)
}

# private function to check the params to mfa
check_params <-function(data, sets) {
  if (!is.data.frame(data) & !is.matrix(data)) {
    stop("data must be a data frame or matrix")
  }
  for (i in 1:length(sets)) {
    if (!is.numeric(sets[[i]]) & !is.character(sets[[i]])) {
      stop("sets must be a numeric or character vector")
    }
  }
  if (length(sets) > ncol(df) | length(sets) < 1) {
    stop("length of sets must be at least 1 and <= the length of data columns")
  }
  return(TRUE)
}

# private function to adjust set values to reduced matrix dimensions
adjust_sets <- function(sets, sets.len) {
  start <- cumsum(sets.len) - sets.len + 1
  end <- cumsum(sets.len)
  sets <- lapply(1:length(sets.len), function(i) start[i]:end[i])
  return(sets)
}

# private function to compute the alpha weight vector
get_alphas <- function(X, sets) {
  # compute alpha value for a submatrix
  alpha_helper <- function(X) {
    X.svd <- svd(X)
    a <- 1 / X.svd$d[1]^2
    return(a)
  }
  # call helper function for each submatrix
  alphas <- sapply(sets, function(i) alpha_helper(X[,i]))
}

# private function to convert the alpha weight vector to a matrix
get_A <- function(a, sets.len) {
  # repeat "a" for the length of each set
  a.rep <- lapply(1:length(a), function(i) rep(a[i], sets.len[i]))
  a.rep <- unlist(a.rep)
  # return as a diagonal matrix
  return(diag(a.rep))
}

# private function to compute the partial factor score array
get_P <- function(X,sets, a, V.tilde) {
  # create an empty array to hold the partial factor scores
  P <- list()
  # fill the array using the provided formula
  for (i in 1:length(sets)) {
    P[[i]] <- 10 * a[i] * X[,sets[[i]]] %*% V.tilde[sets[[i]],]
  }
  return(P)
}

#' @export
print.mfa <- function(x, ...) {
  cat('object "mfa"\n')
  cat(sprintf('eigenvalues: [%d x 1] %.3f %.3f %.3f %.3f %.3f ...',
              length(x$eigenvalues),
              x$eigenvalues[1], x$eigenvalues[2], x$eigenvalues[3], x$eigenvalues[4], x$eigenvalues[5]), "\n")
  cat(sprintf('factor scores: [%d x %d] %.3f %.3f %.3f %.3f %.3f ...',
              nrow(x$factor_scores), ncol(x$factor_scores),
              x$factor_scores[1], x$factor_scores[2], x$factor_scores[3], x$factor_scores[4], x$factor_scores[5]), "\n")
  cat(sprintf('partial factor scores: [%d x 1] %.3f %.3f %.3f %.3f %.3f ...',
              length(x$partial_factor_scores),
              x$partial_factor_scores[[1]][1,1], x$partial_factor_scores[[1]][2,1],
              x$partial_factor_scores[[1]][3,1], x$partial_factor_scores[[1]][4,1],
              x$partial_factor_scores[[1]][5,1]), "\n")
  cat(sprintf('factor loadings: [%d x %d] %.3f %.3f %.3f %.3f %.3f...',
              nrow(x$loadings), ncol(x$loadings),
              x$loadings[1,1], x$loadings[2,1], x$loadings[3,1],
              x$loadings[4,1], x$loadings[5,1]), "\n")
}
