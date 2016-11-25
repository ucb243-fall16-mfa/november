# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}

setwd("../Documents/stat243-fall-2016/problem-sets/final-project/")
df <- read.csv("data/wines.csv", stringsAsFactors = FALSE)

mfa <- function(data, sets=NULL, ncomps = NULL, center = TRUE, scale = TRUE) {
  check_params(data, sets)

  X <- as.matrix(df[,2:54])
  X <- apply(X, 2, function(x) x - mean(x))
  X <- apply(X, 2, function(x) x / sqrt(sum(x * x)))

  start <- c(1,7,13,19,23,30,35,39,45,50)
  end <- c(6,12,18,22,29,34,38,44,49,53)
  loc <- cbind(start,end)
  loc <- loc[sets,]

  a <- get_alphas(X,loc)
  A <- get_A(a,X,loc)
  M <- diag(rep(1 / nrow(X), nrow(X)))

  A.tilde <- sqrt(M) %*% X %*% sqrt(A)
  A.tilde.svd <- svd(A.tilde)
  U.tilde <- solve(sqrt(M)) %*% A.tilde.svd$u
  V.tilde <- solve(sqrt(A)) %*% A.tilde.svd$v
  D.tilde <- A.tilde.svd$d
  loadings <- V.tilde
  factor_scores <- U.tilde %*% diag(D.tilde)

  P <- array(0, c(12,12,10))
  for (i in 1:nrow(loc)) {
    P[,,i] <- 10 * a[i] * X[,loc[i,1]:loc[i,2]] %*% V.tilde[loc[i,1]:loc[i,2],]
  }

  object <- list(
    eigenvalues = D.tilde^2,
    factor_scores = factor_scores,
    partial_factor_scores = P,
    loadings = loadings,
    A = A
  )
  class(object) <- "mfa"
  return(object)
}

# private function to check the params to mfa
check_params <-function(data, sets) {
  if (!is.data.frame(data) & !is.matrix(data)) {
    stop("data must be a data frame or matrix")
  }
  if (!is.numeric(sets) & !is.character(sets)) {
    stop("sets must be a numeric or character vector")
  }
  if (length(sets) > ncol(df) | length(sets) < 1) {
    stop("length of sets must be at least 1 and <= the length of data columns")
  }
  return(TRUE)
}

get_alphas <- function(X,loc) {
  alpha_helper <- function(X) {
    X.svd <- svd(X)
    a <- 1 / X.svd$d[1]^2
    return(a)
  }
  alphas <- apply(loc,1,function(x) alpha_helper(X[,x[1]:x[2]]))
}

get_A <- function(a, X, loc) {
  l <- rep(0,ncol(X))
  for (i in 1:nrow(loc)) {
    s <- loc[i,2] - loc[i,1] + 1
    l[loc[i,1]:loc[i,2]] <- rep(a[i], s)
  }
  return(diag(l))
}

print.mfa <- function(x, ...) {
  cat('object "mfa"\n')
  cat(sprintf('eigenvalues: [%d x 1] %.3f %.3f %.3f %.3f %.3f ...',
              length(x$eigenvalues),
              x$eigenvalues[1], x$eigenvalues[2], x$eigenvalues[3], x$eigenvalues[4], x$eigenvalues[5]), "\n")
  cat(sprintf('factor scores: [%d x %d] %.3f %.3f %.3f %.3f %.3f ...',
              nrow(x$factor_scores), ncol(x$factor_scores),
              x$eigenvalues[1], x$eigenvalues[2], x$eigenvalues[3], x$eigenvalues[4], x$eigenvalues[5]), "\n")
  cat(sprintf('partial factor scores: [%d x %d x %d] %.3f %.3f %.3f %.3f %.3f ...',
              dim(x$partial_factor_scores)[1],
              dim(x$partial_factor_scores)[2],
              dim(x$partial_factor_scores)[3],
              x$partial_factor_scores[1,1,1], x$partial_factor_scores[2,1,1],
              x$partial_factor_scores[3,1,1], x$partial_factor_scores[4,1,1],
              x$partial_factor_scores[5,1,1]), "\n")
  cat(sprintf('factor loadings: [%d x %d] %.3f %.3f %.3f %.3f %.3f...',
              nrow(x$loadings), ncol(x$loadings),
              x$loadings[1,1], x$loadings[2,1], x$loadings[3,1],
              x$loadings[4,1], x$loadings[5,1]), "\n")
}


# plot factor scores
plot.mfa <- function(x, d = c(1, 2), ...) {
  plot(x$factor_scores[,1],x$factor_scores[,2])
  plot(thing$factor_scores[,1],thing$factor_scores[,2], pch=17,
       xlab="", ylab="", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))
  text(thing$factor_scores[,1],thing$factor_scores[,2], labels=df[,1], cex= 0.7, pos=3)
  abline(v=0,lty=2,col="blue")
  abline(h=0,lty=2,col="blue")
}

# plot partial factor scores
plot.mfa <- function(x, d = c(1, 2), ...) {
  par(mfrow=c(2,5))
  for (i in 1:dim(x$partial_factor_scores)[3]) {
    plot(x$partial_factor_scores[,1,i],x$partial_factor_scores[,2,i], pch=17,
         xlab="", ylab="", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), main=paste("Assessor",i))
    text(x$partial_factor_scores[,1,i],x$partial_factor_scores[,2,i], labels=df[,1], cex= 0.7, pos=3)
    abline(v=0,lty=2,col="blue")
    abline(h=0,lty=2,col="blue")
  }
}

extract <- function(mfa) {
  sing_vals <- round(sqrt(mfa$eigenvalues),3)
  eigen <- round(mfa$eigenvalues,3)
  eigen_cum <- round(cumsum(mfa$eigenvalues),3)
  perc_inertia <- round(mfa$eigenvalues / sum(mfa$eigenvalues),2) * 100
  perc_cum <- round(cumsum(mfa$eigenvalues / sum(mfa$eigenvalues)),2) * 100
  df <- cbind(sing_vals,eigen,eigen_cum,perc_inertia,perc_cum)
}

obs_to_dim <- function(mfa) {
  M <- sapply(1:length(mfa$eigenvalues), function(i) (1/12 * mfa$factor_scores[,i]^2) / mfa$eigenvalues[i])
  return(M)
}

var_to_dim <- function(mfa) {
  mfa$A %*% mfa$loadings^2
}

tbl_to_dim <- function(mfa,loc) {
  v2d <- var_to_dim(mfa)
  M <- matrix(nrow=10,ncol=12)
  for (i in 1:nrow(loc)) {
    M[i,] <- apply(v2d,2,function(x) sum(x[loc[i,1]:loc[i,2]]))
  }
  return(M)
}

RV <- function(table1, table2) {
  num <- sum(diag((tcrossprod(table1) %*% tcrossprod(table2))))
  den <- sqrt(sum(diag(tcrossprod(table1) %*% tcrossprod(table1))) *
              sum(diag(tcrossprod(table2) %*% tcrossprod(table2))))
  return(num/den)
}

RV_table <- function(dataset, sets=list(1:3,4:5,6:10)) {
  M <- diag(length(sets))
  for (i in 1:length(sets)) {

  }

}
