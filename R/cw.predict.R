cw.predict <-
function(Xnew, res.cw){
 

  # -----------------------------------------------------------------------
  #  1. Required data and parameters
  # -----------------------------------------------------------------------
  
  Xnew    <- as.matrix(Xnew)                   
  appel   <- as.list(res.cw$call) 
  X.train <- as.matrix(eval.parent(appel$X))   
  Y.train <- as.matrix(eval.parent(appel$Y))
  G       <- as.matrix(eval.parent(appel$G))
  P       <- dim(Xnew)[2]
  Q       <- dim(Y.train)[2]
  
  
  # ---------------------------------------------------------------------------
  # 1bis. Preliminary tests
  # ---------------------------------------------------------------------------
  
  if (any(is.na(Xnew))) 
    stop("No NA values are allowed")
  if (any(colnames(Xnew) != colnames(X.train)))
    stop("Xnew and X must be measured on the same variables")
  if (!inherits(res.cw, "cwmultiblock")) 
    stop("Object of type 'cwmultiblock' expected")
  
  
  #------------------------------------------------------------------------------
  # 2. Pre-processing of Xnew and X.train (depending on X train)
  #------------------------------------------------------------------------------
  
  Xnew.c     <- sweep(Xnew, 2, colMeans(X.train), FUN="-")
  Xnew.cr    <- sweep(Xnew.c, 2, apply(X.train, 2, sd), FUN="/")
  X.train.cr <- scale(X.train, center = TRUE, scale = TRUE)
  
  
  # -------------------------------------------------------------------
  # 3. No expected cluster structure (G = 1)
  # -------------------------------------------------------------------
  
  if (G == 1){
    
    res <- list()
    
    # 21. Prediction of the dependent variable values according to the normalized train data
    res$Ypred.cr            <- sapply(1:Q, function(q) Xnew.cr %*% res.cw$beta.cr[2:(P+1), q])
    colnames(res$Ypred.cr)  <- colnames(Y.train)
    rownames(res$Ypred.cr)  <- rownames(Xnew.cr)
    
    # 22. Prediction of the dependent variable values according to the raw train data
    res$Ypred.raw           <- sapply(1:Q, function(q) matrix(rep(res.cw$beta.raw[1, q], each=dim(Xnew)[1], nrow=dim(Xnew)[1])) + Xnew %*% res.cw$beta.raw[2:(P+1), q])
    colnames(res$Ypred.raw) <- colnames(Y.train)
    rownames(res$Ypred.raw) <- rownames(Xnew)
    

    
  # -------------------------------------------------------------------
  # 4. Expected cluster structure (G > 1)
  # ------------------------------------------------------------------- 
    
  } else {
    
    # 41. Affectation of new observations (KNN rule)
    ZX.train           <- as.data.frame(cbind(res.cw$cluster, X.train.cr))
    ZX.train[, 1]      <- as.factor(ZX.train[, 1])
    ZXnew              <- as.data.frame(cbind(rep(1, dim(Xnew.cr)[1]), Xnew.cr))
    ZXnew[, 1]         <- as.factor(ZXnew[, 1])
    names(ZX.train)[1] <- names(ZXnew)[1] <- c("Z")
    fit.train          <- train.kknn(Z ~ ., ZX.train, kmax = min(25, (dim(ZX.train)[1]/2)), kernel = c("rectangular", "triangular", "epanechnikov", "gaussian", "rank", "optimal"))
    # fit.train          <- train.kknn(Z ~ ., ZX.train, kmax = 8, kernel = c("rectangular", "triangular", "epanechnikov", "gaussian", "rank", "optimal"))
    res.kknn           <- kknn(Z ~., train = as.data.frame(ZX.train), test = as.data.frame(ZXnew), 
                         k = fit.train$best.parameters$k, kernel = fit.train$best.parameters$kernel)
    clusternew         <- res.kknn$fitted.values
    names(clusternew)  <- row.names(Xnew)
    G.new              <- as.numeric(levels(factor(clusternew)))
    
    
    # 42. Prediction of the dependent variable values according to the normalized train data
    Xnew.cr.g <- lapply(1:G, function(g)  Xnew.cr[which(clusternew == g), , drop = FALSE])
    Ypred.cr  <- list()
    for (g in G.new){
      Ypred.cr[[g]]           <- sapply(1:Q, function(q) matrix(rep(res.cw$beta.cr[[g]][1, q], each=dim(Xnew.cr.g[[g]])[1]), nrow=dim(Xnew.cr.g[[g]])[1]) + Xnew.cr.g[[g]] %*% res.cw$beta.cr[[g]][2:(P+1), q])
      if (is.vector(Ypred.cr[[g]])){Ypred.cr[[g]] <- t(as.matrix(Ypred.cr[[g]]))}
      colnames(Ypred.cr[[g]]) <- colnames(Y.train)
      rownames(Ypred.cr[[g]]) <- rownames(Xnew.cr.g[[g]])
    }
    
    
    # 43. Prediction of the dependent variable values according to the raw train data
    Xnew.g    <- lapply(1:G, function(g) Xnew[which(clusternew == g), , drop = FALSE])
    Ypred.raw <- list()
    for (g in G.new){
      Ypred.raw[[g]]           <- sapply(1:Q, function(q) matrix(rep(res.cw$beta.raw[[g]][1, q], each=dim(Xnew.g[[g]])[1]), nrow=dim(Xnew.g[[g]])[1]) + as.matrix(Xnew.g[[g]]) %*% res.cw$beta.raw[[g]][2:(P+1), q])
      if (is.vector(Ypred.raw[[g]])){Ypred.raw[[g]] <- t(as.matrix(Ypred.raw[[g]]))}
      colnames(Ypred.raw[[g]]) <- colnames(Y.train)
      rownames(Ypred.raw[[g]]) <- rownames(Xnew.g[[g]])
    }
    
    # 44. Result storage
    res             <- list()
    res$clusternew  <- clusternew
    res$Ypred.cr    <- Ypred.cr
    res$Ypred.raw   <- Ypred.raw
  }
  
  return(res)
}
