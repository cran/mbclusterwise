cw.multiblock <-
function(Y, X, blo, option=c("none", "uniform"), G, H, INIT = 20, method = c("mbpls", "mbpcaiv", "mbregular"), Gamma = NULL, parallel.level = c("high", "low")){
  
  
  # ---------------------------------------------------------------------------
  # 0. Preliminary tests
  # ---------------------------------------------------------------------------
  
  if (any(is.na(Y)) | any(is.na(X))) 
    stop("No NA values are allowed")
  if (is.null(rownames(X)) | is.null(rownames(Y))) 
    stop("X and Y must have rownames")
  if (is.null(colnames(X))) 
    colnames(X) <- paste("X", 1:dim(X)[2], sep=".")
  if (is.null(colnames(as.data.frame(Y)))) 
    colnames(Y) <- paste("Y", 1:dim(as.data.frame(Y))[2], sep=".")
  if (any(row.names(X) != row.names(Y)))
    stop("X and Y must have the same rows")
  if (sum(blo) != ncol(X)) 
    stop("Non convenient blocks parameter")
  if (nrow(X) < 5)
    stop("Minimum five observations are required")
  if (min(blo) < 2)
    stop("Minimum two variables per explanatory block are required")
  if (INIT < 2)
    stop("Minimum two initializations are required")
  if (method == "mbpls" | method == "mbpcaiv"){
    if (!is.null(Gamma)) stop("Non convenient value for parameter Gamma")
  } else {
    if (is.null(Gamma)) stop("Non convenient value for parameter Gamma")
  }
  option <- match.arg(option)
  method <- match.arg(method)
  
  
  # -----------------------------------------------------------------------
  #  1. Required parameters
  # -----------------------------------------------------------------------
  
  N <- nrow(X)
  P <- ncol(X)
  Q <- ncol(as.data.frame(Y))
  K <- length(blo)
  
  
  #------------------------------------------------------------------------------
  # 2. Preliminary function
  #------------------------------------------------------------------------------
  

  mbreg <- function(y, x, h, gamma){
      dudiy <- dudi.pca(y, center = FALSE, scale = FALSE, scannf = FALSE)
      ktabx <- ktab.data.frame(df=data.frame(x), blocks = blo, tabnames = paste("Tab", c(1:K), sep="."))
      if (method == "mbpls"){
        gamma <- NULL
        m <- mbpls.fast(dudiy, ktabx, scale = FALSE, option = option, H = h)
      } else if (method == "mbpcaiv"){
        gamma <- NULL
        m <- mbpcaiv.fast(dudiy, ktabx, scale = FALSE, option = option, H = h)
      } else {
        m <- mbregular(dudiy, ktabx, scale = FALSE, option = option, H = h, gamma = Gamma)
      }
      h.opt <- min(h, dim(as.data.frame(m$lX))[2])
      for (l in 1:h.opt){
        XYcoef    <- sapply(m$XYcoef, function(x) x[, l])
        intercept <- sapply(m$intercept, function(x) x[, l])
        fitted    <- sapply(m$fitted, function(x) x[, l])
      }
      return(list(co = rbind(intercept, XYcoef), err = m$crit.reg[h.opt], pred = fitted, hopt = h.opt))  
    }
  
  
  
  #------------------------------------------------------------------------------
  # 3. Pre-processing of X and Y
  #------------------------------------------------------------------------------
  
  # Raw variable means and variances
  meanY  <- apply(as.matrix(Y), 2, mean)
  sdY    <- apply(as.matrix(Y), 2, sd)
  meanX  <- apply(as.matrix(X), 2, mean)
  sdX    <- apply(as.matrix(X), 2, sd)
  X.raw  <- as.matrix(X)
  
  # Data normalization
  Y <- scale(Y, center = TRUE, scale = TRUE)
  X <- scale(X, center = TRUE, scale = TRUE)
  
  
  #------------------------------------------------------------------------------
  # 4. Multiblock analysis (when no expected cluster ; G=1)
  #------------------------------------------------------------------------------
  
  if (G == 1){
    res.mbreg               <- mbreg(Y, X, H, gamma = Gamma)
    res                     <- list()
    res$error               <- res.mbreg$err   # NB : performed on CR data
    res$beta.cr             <- res.mbreg$co
    res$beta.raw            <- matrix(NA, nrow = P+1, ncol = Q, dimnames = list(c("intercept", colnames(X)), colnames(Y)))
    res$beta.raw[-1, ]      <- res.mbreg$co[-1, , drop = FALSE] * matrix(rep(sdY, each=P), ncol=Q) / t(matrix(rep(sdX, each=Q), nrow=Q))
    res$beta.raw[1, ]       <- meanY - meanX %*% res$beta.raw[-1, ]
    res$hopt                <- res.mbreg$h.opt
    res$Ypred.cr            <- res.mbreg$pred
    res$Ypred.raw           <- sapply(1:Q, function(q) matrix(rep(res$beta.raw[1, q], each=N, nrow=N)) + X.raw %*% res$beta.raw[2:(P+1), q])
    colnames(res$Ypred.raw) <- colnames(Y)
    rownames(res$Ypred.raw) <- rownames(X.raw)
    
    
    
    
  #---------------------------------------------------------------------------------------
  # 5. Sequential clusterwise multiblock analysis (when expected cluster structure ; G>1)
  #---------------------------------------------------------------------------------------
  } else {
    
    # Preparation of the parallelized processing
    nodes <- detectCores()
    if (parallel.level == "low"){nodes <- 2}
    cl    <- makeCluster(nodes)
    set.seed((as.numeric(Sys.time()) - floor(as.numeric(Sys.time()))) * 1e8 -> seed)
    clusterSetRNGStream(cl, iseed = seed)
    registerDoParallel(cl)
    
    
    # Preparation of the result storage
    crit.stock      <- list()
    cluster.stock   <- list()
    init            <- NULL
    
    res.init = foreach(init = 1:INIT, .export=c('mbreg', 'mbpls.fast', 'mbpcaiv.fast', 'mbregular'), .packages=c('ade4')) %dopar%{
      
      to.res  <- NULL
      crit    <- rep(0, N)
      
      #------------------------------------------------------------------------------
      # 51. Several initializations of observations into G clusters
      #------------------------------------------------------------------------------
      
      # Four optimal initializations and then random ones
      if (init == 1){
        cluster <- kmeans(X, G, iter.max = 30, nstart = 10, algorithm = "Forgy")$cluster
        if (min(table(cluster)) <= 2){cluster <- 1 + floor(runif(N)*G)}
      } else if (init == 2){
        cluster <- kmeans(cbind(Y, X), G, iter.max = 30, nstart = 10, algorithm = "Forgy")$cluster
        if (min(table(cluster)) <= 2){cluster <- 1 + floor(runif(N)*G)}
      } else if (init == 3){
        dudiy   <- dudi.pca(Y, center = FALSE, scale = FALSE, scannf = FALSE)
        ktabx   <- ktab.data.frame(df = data.frame(X), blocks = blo, tabnames = paste("Tab", c(1:K), sep="."))
        pls.tot <- mbpls.fast(dudiy, ktabx, scale = FALSE, option = option, H = H)
        cluster <- kmeans(as.data.frame(pls.tot$lX)[, 1:min(H, dim(as.data.frame(pls.tot$lX))[2])], G, iter.max = 30, nstart = 10, algorithm = "Forgy")$cluster
        if (min(table(cluster)) <= 2){cluster <- 1 + floor(runif(N)*G)}
      } else if (init == 4){
        dudiy       <- dudi.pca(Y, center = FALSE, scale = FALSE, scannf = FALSE)
        ktabx       <- ktab.data.frame(df = data.frame(X), blocks = blo, tabnames = paste("Tab", c(1:K), sep="."))
        mbpcaiv.tot <- mbpcaiv.fast(dudiy, ktabx, scale = FALSE, option = option, H = H)
        cluster     <- kmeans(as.data.frame(mbpcaiv.tot$lX)[,  1:min(H, dim(as.data.frame(mbpcaiv.tot$lX))[2])], G, iter.max = 30, nstart = 10, algorithm = "Forgy")$cluster
        if (min(table(cluster)) <= 2){cluster <- 1 + floor(runif(N)*G)}
      } else {
        cluster <- 1 + floor(runif(N)*G) 
      }
      
      # Avoid too small clusters
      while (min(table(cluster)) <= 2){cluster <- 1 + floor(runif(N)*G)}
      
      
      
      #------------------------------------------------------------------------------
      # 52. Sequential assignation of each randomly selected observation
      #------------------------------------------------------------------------------
      
      N.random <- sample(c(1:N), N, replace = FALSE)
      
      for(n in 1:N){   # Loop for randomly selected observations
        
        n.sel <- N.random[n]  
        
        # -------------------------------------------------------------------------------------------
        # a. Compute G separate multiblock analyses for the selected observation put in each cluster
        #--------------------------------------------------------------------------------------------
        
        mat.error    <- matrix(NA, nrow = G, ncol = G)
        error.decide <- rep(NA, times = G)
        cluster.n    <- cluster
        for(g in 1 : G){
          # Solution for n which does not belong to cluster g
          vec.mg           <- rep(1:G)[-g]
          cluster.n[n.sel] <- vec.mg[1]
          index.n          <- which(cluster.n == g)
          if (length(index.n) == 1){break}
          local.n         <- mbreg(as.data.frame(Y)[index.n, ], X[index.n, ], H, Gamma)
          mat.error[, g]  <- local.n$err

          # Solution for n which belongs to cluster g
          cluster.n[n.sel] <- g
          index.n          <- which(cluster.n == g)
          if (length(index.n) == 1){break}
          local.n         <- mbreg(as.data.frame(Y)[index.n, ], X[index.n, ], H, Gamma)
          mat.error[g, g] <- local.n$err
          
          error.decide <- rowSums(mat.error)
        }
        
        
        # --------------------------------------------------------------------------------------
        # b. New assignation of the selected observation
        #--------------------------------------------------------------------------------------
        
        if (any(is.na(error.decide)) == FALSE){
          cluster[n.sel] <- which.min(error.decide)
          crit[n]        <- min(error.decide, na.rm = TRUE)
        }
        
      } # End of the loop for observations
      
      
      # Decreasing criterion and associated clustering for each intialization
      to.res$crit.stock    <- crit
      to.res$cluster.stock <- cluster
      to.res
      
    } # End of the loop for initializations
    stopCluster(cl)
    
    
    #---------------------------------------------------------------------------------------
    # 6. Optimal solution
    #---------------------------------------------------------------------------------------
    
    # Selection of the best intialization
    crit.stock           <- lapply(1:INIT, function(x) res.init[[x]]$crit.stock)
    cluster.stock        <- lapply(1:INIT, function(x) res.init[[x]]$cluster.stock)
    init.sel             <- which.min(unlist(lapply(crit.stock, min)))
    crit.final           <- crit.stock[[init.sel]]
    cluster.final        <- cluster.stock[[init.sel]]  
    names(cluster.final) <- rownames(X)
    
    
    # Processing of the clusterwise multiblock models associated with the best initialization
    local <- lapply(1:G, function(g) mbreg(as.data.frame(Y)[which(cluster.final == g), ], X[which(cluster.final == g), ], H, gamma = Gamma))
    beta  <- lapply(1:G, function(g) local[[g]]$co)

    
    # Correction of the regression coefficients and the predicted Y values according to raw data
    X.raw.g  <- lapply(1:G, function(g) X.raw[which(cluster.final == g), , drop = FALSE])
    beta.raw <- lapply(1:G, function(g) matrix(NA, nrow=P+1, ncol=Q, dimnames = list(c("intercept", colnames(X)), colnames(Y))))
    pred.raw <- list()
    for(g in 1: G){
      beta.raw[[g]][-1, ]     <- beta[[g]][-1, ] * matrix(rep(sdY, each=P), ncol=Q) / t(matrix(rep(sdX, each=Q), nrow=Q))
      beta.raw[[g]][1, ]      <- beta[[g]][1, ] * sdY + meanY - meanX %*% beta.raw[[g]][-1, ]
      pred.raw[[g]]           <- sapply(1:Q, function(q) matrix(rep(beta.raw[[g]][1, q], each=dim(X.raw.g[[g]])[1]), nrow=dim(X.raw.g[[g]])[1]) + X.raw.g[[g]] %*% beta.raw[[g]][2:(P+1), q])
      colnames(pred.raw[[g]]) <- colnames(Y)
      rownames(pred.raw[[g]]) <- rownames(X.raw.g[[g]])
    }
    
    
    #------------------------------------------------------------------------------
    # 7. Result storage
    #------------------------------------------------------------------------------    
    
    res           <- list()
    res$error     <- crit.final   # NB : performed on CR data
    res$beta.cr   <- beta
    res$beta.raw  <- beta.raw
    res$cluster   <- cluster.final
    res$hopt      <- lapply(1:G, function(g) local[[g]]$hopt)
    res$Ypred.cr  <- lapply(1:G, function(g) local[[g]]$pred)
    res$Ypred.raw <- pred.raw
  }
  
  res$call   <- match.call()
  class(res) <- c("cwmultiblock")
  return(res)
  
}
