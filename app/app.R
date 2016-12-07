#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Multiple Factor Analysis"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("checkGroup", label = h3("Assessor"),
                         choices = list("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5,
                                        "6" = 6, "7" = 7, "8" = 8, "9" = 9, "10" = 10),
                         selected = c(1:10)),
      radioButtons("radio", label = h3("Plot Type"),
                   choices = list("compromise" = 1, "partial factor scores" = 2, "variable loadings" = 3),
                   selected = 1)
        ),
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
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

  mfa <- function(df, sets, ncomps, center = TRUE, scale = TRUE) {
    # call helper function to check that params are in proper format
    # get descriptor vector

    names <- names(df)[unlist(sets)]
    # pull out data columns specified in sets
    X <- sapply(sets, function(x) df[,x][])

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


  plot_comp <- function(mfa, d = c(1, 2), ...) {
    par(mfrow=c(1,1))
    plot(mfa$factor_scores[,1],mfa$factor_scores[,2], pch=17,
         xlab="", ylab="", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), main="compromise of the tables")
    text(mfa$factor_scores[,1],mfa$factor_scores[,2], labels=mfa$labels, cex= 0.7, pos=3)
    abline(v=0,lty=2,col="blue")
    abline(h=0,lty=2,col="blue")
  }

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
  output$distPlot <- renderPlot({
    if (!is.null(input$checkGroup)) {
      sets=list(2:7,8:13,14:19,20:23,24:30,31:35,36:39,40:45,46:50,51:54)
      selected <- sets[as.numeric(input$checkGroup)]
      df <- read.csv("wines.csv", stringsAsFactors = FALSE)
      mfa1 <- mfa(df, sets=selected)
      if (input$radio == 1) plot_comp(mfa1)
      if (input$radio == 2) plot_pfs(mfa1)
      if (input$radio == 3) plot_loadings(mfa1)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)

