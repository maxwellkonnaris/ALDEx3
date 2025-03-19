##' ALDEx3 Linear Models
##'
##' 
##' @title ALDEx3 Linear Modules
##' @param Y An (D x N) matrix of sequence counts where D is the number of features (taxa/genes) and N is the number of samples
##' @param X Either a formula (requires non-null data parameter) or a model matrix (P x N) where P is the number of linear model covariates
##' @param data Data frame containing variables in formula X (must have N rows)
##' @param nsample Number of Monte Carlo replicates for Dirichlet sampling
##' @param GAMMA Scale model specification. Can be either:
##' \itemize{
##'   \item A function taking arguments (X, Y, logWpara) that returns an (N x nsample) matrix of scale factors on log2 scale
##'   \item An (N x nsample) matrix of pre-computed scale factors
##' }
##' @param streamsize Approximate memory footprint in Mb when using streaming (default 8000). Set to Inf to disable streaming.
##'   This should be thought of as the number of Mb for each
##'   streaming chunk. If D*N*nsample*8/1000000 is less than streamsize then no
##'   streaming will be performed. Note, to conserve memory, samples from the
##'   Dirichlet and scale models will not be returned if streaming is used.
##' @param return.samples (default TRUE) if true, return samples for logWpara
##'   composition and logWperp (scale). Will override to FALSE if streaming is
##'   required.
##' @param p.adjust.method Multiple testing correction method (default "BH"). See \code{p.adjust} for options.
##' @return A list containing:
##' \itemize{
##'   \item estimate - (P x D x nsample) array of coefficient estimates
##'   \item std.error - (P x D x nsample) array of standard errors
##'   \item p.val - (P x D) matrix of adjusted p-values
##' }
##' TODO p.value calcluation may be slightly different than in current ALDEx3 -- need to check.
##' 
##' @export
##' @examples
##' set.seed(43254)
##' 
##' # Simulation parameters
##' N <- 300
##' D <- 100
##' DE <- 40
##' mc.samples <- 200
##' 
##' # Create metadata (2-group design)
##' metadata <- data.frame(condition = rep(c("A", "B"), each = N/2))
##' 
##' # Simulate baseline abundances
##' sim_A <- matrix(rnorm(D*N, mean = 5, sd = 2), nrow = D, ncol = N)
##' 
##' # Add differential abundance effects
##' DE_taxa <- sample(1:D, DE)
##' lfcs <- rnorm(DE, mean = 1, sd = 0.5)
##' design_matrix <- model.matrix(~condition, metadata)
##' sim_A[DE_taxa,] <- sim_A[DE_taxa,] + lfcs %*% t(design_matrix[,2])
##' 
##' # Generate count data
##' sim_Y <- apply(2^sim_A, 2, function(x) rmultinom(1, 1e6, prob = x))
##' 
##' # Define custom scale model function
##' gamma_func <- function(X, Y, logWpara) {
##'   # Create N x nsample matrix of scale factors
##'   matrix(rnorm(N*mc.samples, mean = 0.5, sd = 0.5), nrow = N, ncol = mc.samples)
##' }
##' 
##' # Run ALDEx3 analysis
##' aldex3.res <- aldex.lm(sim_Y, ~condition, data = metadata, 
##'                       nsample = mc.samples, GAMMA = gamma_func)
##'                       
##' # Extract mean coefficients
##' coef_means <- apply(aldex3.res$estimate, c(1,2), mean)
##' @export
##' @author Justin Silverman
aldex.lm <- function(Y, X, data=NULL, nsample=2000,  GAMMA=NULL,
                     streamsize=8000, return.samples=FALSE,
                     p.adjust.method="BH") {
  N <- ncol(Y)
  D <- nrow(Y)

  ## calculate streaming threshold
  stream <- FALSE
  if (N*D*nsample*8/1000000 > streamsize) stream <- TRUE

  ## compute model matrix
  if (inherits(X, "formula")) {
    if (is.null(data)) stop("data should not be null if X is a formula")
    
    X <- t(model.matrix(X, data))
    
    ## Error checking
    if (!is.matrix(X)) stop("model.matrix() did not return a matrix")
    if (ncol(X) != ncol(Y)) stop(paste("Mismatch: X must have", ncol(Y), "rows to match Y's samples, but has", ncol(X)))
    if (nrow(X) == 0) stop("Error: model matrix X has zero columns, ensure your formula includes valid predictors.")
  } else {
    ## Ensure X is a matrix
    if (!is.matrix(X)) stop("X must be a formula or a design matrix")
    
    ## Check dimensions
    if (ncol(X) != ncol(Y)) stop(paste("Mismatch: X must have", ncol(Y), "rows to match Y's samples, but has", ncol(X)))
    if (nrow(X) == 0) stop("Error: Design matrix X has zero columns, ensure it contains valid predictors.")
  }
  
  ## Ensure numeric values in X
  if (!all(is.numeric(X))) stop("Error: X contains non-numeric values. Ensure all covariates are numerical.")
  
  ## perform streaming 
  out <- list()
  nsample.remaining <- nsample
  iter <- 1
  if (stream) {
    nsample.local <- floor(streamsize*1000000/(N*D*8))
    if (nsample.local < 1) stop("streamsize too small")
  } else {
    nsample.local <- nsample
  }
  while (nsample.remaining > 0) {
    nsample.remaining <- nsample.remaining - nsample.local
    out[[iter]] <- aldex.lm.internal(Y, X, nsample.local, GAMMA, stream)
    iter <- iter+1
  }
  ## combine output of the different streams
  out <- combine.streams(out)

  # p-value calculations, accounting for sign changes
  p.lower <- apply(2*out$p.lower, c(1,2,3), function(item) min(1, item))
  p.upper <- apply(2*out$p.upper, c(1,2,3), function(item) min(1, item))
  p.lower.adj <- apply(p.lower, c(1,3), function(item) {
    p.adjust(item, method=p.adjust.method)
  })
  p.upper.adj <- apply(p.upper, c(1,3), function(item) {
    p.adjust(item, method=p.adjust.method)
  })
  p.lower.mean <- apply(p.lower, c(2,1), mean)
  p.upper.mean <- apply(p.upper, c(2,1), mean)
  p.res <- c()
  for(col_i in 1:ncol(p.lower.mean)) {
    tmp_mat <- cbind(p.lower.mean[,col_i],
                     p.upper.mean[,col_i])
    p.res <- rbind(p.res, apply(tmp_mat, 1, min))
  }
  p.lower.mean.adj <- apply(p.lower.adj, c(1,2), mean)
  p.upper.mean.adj <- apply(p.upper.adj, c(1,2), mean)
  p.adj.res <- c()
  for(col_i in 1:ncol(p.lower.mean.adj)) {
    tmp_mat <- cbind(p.lower.mean.adj[,col_i],
                     p.upper.mean.adj[,col_i])
    p.adj.res <- rbind(p.adj.res, apply(tmp_mat, 1, min))
  }

  return(list(estimate=out$estimate, std.error=out$std.error, p.val=p.res,
              p.val.adj=p.adj.res))
  ## TODO write a good "summary" function and wrap this all in S3 class
}


aldex.lm.internal <- function(Y, X, nsample, GAMMA=NULL, stream) {
  N <- ncol(Y)
  D <- nrow(Y)
  
  ## dirichlet sample
  logWpara <- log2(rDirichletMat(nsample, Y+0.5))

  ## sample from scale model
  if (is.null(GAMMA)){
    stop("You probably want a scale model :)")
  } else if (is.function(GAMMA)) {
    logWperp <- GAMMA(X, Y, logWpara)
  } else if (is.matrix(GAMMA)) { 
    logWperp <- GAMMA
  }
  
  ## Correct Dimensions Check
  expected_dim <- c(ncol(Y), nsample)  
  actual_dim <- dim(logWperp)
  
  if (!all(actual_dim == expected_dim)) {
    stop(paste0("GAMMA function returned wrong dimensions. Expected: ", 
                paste(expected_dim, collapse=" x "), 
                " received: ", paste(actual_dim, collapse=" x ")))
  }
  
  # Reasonable biological range check for log2 
  if (any(abs(logWperp) > 10)) {
    warning(paste(
      "|log2 scale factors| >10 detected in", sum(abs(logWperp) > 10), "cases.",
      "Verify these are intentional and that GAMMA function outputs 
      or GAMMA matrix contains log2 values."
    ))
  }
  
  # Natural log mistake check (ln(2) ≈ 0.693 threshold)
  if (median(abs(logWperp)) < 1 && max(abs(logWperp)) > 5) {
    warning(paste(
      "Suspicious scale distribution detected.",
      "Are you certain GAMMA values are log2 (not natural log)?",
      "Median:|", round(median(abs(logWperp)), 2), "| Max:|", round(max(abs(logWperp)), 2), "|"
    ))
  }

  ## compute scaled abundances (W)
  logW <- sweep(logWpara, c(2,3), logWperp, FUN=`+`)

  ## fit linear model
  out <- fflm(aperm(logW, c(2,1,3)),t(X)) # TODO change fflm so it gives correct
                                          # output dimensions and takes correct
                                          # inputs dimensions -- without needing
                                          # aperm
  if (!stream){
    out$logWpara <- logWpara
    out$logWperp <- logWperp
  }
  return(out)
 }
