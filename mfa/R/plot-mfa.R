#' @title plot_comp
#' @description Plots the compromise of an mfa object
#' @param mfa an object of class \code{"mfa"}
#' @param d a vector of dimensions to display
#' @param \dots arguments to be passed to/from other methods
#' @export
#' @examples
#'  \dontrun{
#'  # create an mfa object
#'  mfa1 <- mfa(data, sets=list(2:7,8:13,14:19)
#'
#'  plot_comp(mfa1)
#'  }
plot_comp <- function(mfa, d = c(1, 2), ...) {
  par(mfrow=c(1,1))
  plot(mfa$factor_scores[,1],mfa$factor_scores[,2], pch=17,
       xlab="", ylab="", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), main="compromise of the tables")
  text(mfa$factor_scores[,1],mfa$factor_scores[,2], labels=mfa$labels, cex= 0.7, pos=3)
  abline(v=0,lty=2,col="blue")
  abline(h=0,lty=2,col="blue")
}

#' @title plot_pfs
#' @description Plots the partial factor scores of an mfa object
#' @param mfa an object of class \code{"mfa"}
#' @param d a vector of dimensions to display
#' @param \dots arguments to be passed to/from other methods
#' @export
#' @examples
#'  \dontrun{
#'  # create an mfa object
#'  mfa1 <- mfa(data, sets=list(2:7,8:13,14:19)
#'
#'  plot_pfs(mfa1)
#'  }
# plot partial factor scores
plot_pfs <- function(mfa, d = c(1, 2), ...) {
  mfval <- ifelse(length(mfa$sets) >= 5, par(mfrow=c(2,5)), par(mfrow=c(1,length(mfa$sets))))
  for (i in 1:length(mfa$sets)) {
    # plot partial factor scores
    plot(mfa$partial_factor_scores[[i]][,1],mfa$partial_factor_scores[[i]][,2], pch=17,
         xlab="", ylab="", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), main=paste("Assessor",i))
    text(mfa$partial_factor_scores[[i]][,1],mfa$partial_factor_scores[[i]][,2], labels=mfa$labels, cex= 0.7, pos=3)
    M <-mfa$loadings[1:2,mfa$sets[[i]]]
    par(new = TRUE)
    plot(M[1,], M[2,], pch=22, bg=22, xlab="", ylab="", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))
    text(M[1,], M[2,], labels=mfa$names[mfa$sets[[i]]], cex= 0.7, pos=3)
    abline(v=0,lty=2,col="blue")
    abline(h=0,lty=2,col="blue")
  }
  par(mfrow=c(1,1))
}

#' @title plot_loadings
#' @description Plots the loadings of an mfa object
#' @param mfa an object of class \code{"mfa"}
#' @param d a vector of dimensions to display
#' @param \dots arguments to be passed to/from other methods
#' @export
#' @examples
#'  \dontrun{
#'  # create an mfa object
#'  mfa1 <- mfa(data, sets=list(2:7,8:13,14:19)
#'
#'  plot_loadings(mfa1)
#'  }
# plot loadings
plot_loadings <- function(mfa, d = c(1, 2), ...) {
  par(mfrow=c(2,5))
  for (i in 1:length(mfa$sets)) {
    M <-mfa$loadings[1:2,mfa$sets[[i]]]
    plot(M[1,], M[2,], pch=22, bg=22, xlab="", ylab="", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), main=paste("Assessor",i))
    text(M[1,], M[2,], labels=mfa$names[mfa$sets[[i]]], cex= 0.7, pos=3)
    abline(v=0,lty=2,col="blue")
    abline(h=0,lty=2,col="blue")
  }
  par(mfrow=c(1,1))
}
