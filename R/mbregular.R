mbregular <-
function(dudiY, ktabX, scale=FALSE, option = c("none", "uniform"), H, gamma){

  
  
  # ---------------------------------------------------------------------------
  # 0. Preliminary tests
  # ---------------------------------------------------------------------------
  
  if (!inherits(dudiY, "dudi"))
    stop("object 'dudi' expected")
  if (!inherits(ktabX, "ktab"))
    stop("object 'ktab' expected")
  if (!(is.logical(scale)))
    stop("Non convenient selection for scaling")
  option <- match.arg(option)
  
  
  ## -------------------------------------------------------------------------------
  ##			Arguments and data transformation
  ## -------------------------------------------------------------------------------
  
  ginv <- function(X, tol = sqrt(.Machine$double.eps)){
    if (!is.matrix(X)) 
      X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
      Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) 
      Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
      array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
  }
  
  ## Raw variable means and variances
  meanY        <- colMeans(as.matrix(dudiY$tab))
  sdY          <- apply(as.matrix(dudiY$tab), 2, sd)
  nblo         <- length(ktabX$blo)
  meanX        <- colMeans(cbind.data.frame(lapply(unclass(ktabX)[1 : nblo], scale, center = FALSE, scale = FALSE)))
  names(meanX) <- col.names(ktabX)
  sdX          <- apply(cbind.data.frame(lapply(unclass(ktabX)[1 : nblo], scale, center = FALSE, scale = FALSE)), 2, sd)
  
  ## Preparation of the data frames
  nr    <- nrow(as.matrix(dudiY$tab)) 
  ncolY <- ncol(as.matrix(dudiY$tab)) 
  Y     <- as.matrix(scale(dudiY$tab, center = TRUE, scale = scale))
  Xk    <- lapply(unclass(ktabX)[1 : nblo], scalewt, wt = ktabX$lw, center = TRUE, scale = scale)      # X data with biaised variance
  if (scale == TRUE){Xk <- lapply(1:nblo, function(k) Xk[[k]] * sqrt((nr-1)/nr))}                      # X data with unbiaised variance
  
  ## Block weighting
  if (option[1] == "uniform"){
    In.Y <- sqrt((1/(nr-1)) * sum(diag(crossprod(Y))))
    Y    <- Y / In.Y
    In.Xk <- list()
    for (k in 1 : nblo){
      In.Xk[[k]] <- sqrt((nblo/(nr-1)) * sum(diag(crossprod(Xk[[k]]))))
      Xk[[k]]    <- Xk[[k]] / sqrt((nblo/(nr-1)) * sum(diag(crossprod(Xk[[k]]))))
    }
  }
  X           <- cbind.data.frame(Xk)
  colnames(X) <- col.names(ktabX)
  ncolX       <- ncol(X)
  maxdim      <- H
  
  
  ##-----------------------------------------------------------------------
  ##                         Prepare the outputs
  ##-----------------------------------------------------------------------
  
  dimlab <- paste("Ax", 1:maxdim, sep = "")
  res    <- list(lX     = matrix(0, nrow = nr, ncol = maxdim, dimnames = list(row.names(dudiY$tab), dimlab)), 
                 XYcoef    = lapply(1:ncolY, function(q)  matrix(0, nrow = ncolX, ncol = maxdim, dimnames = list(colnames(X), dimlab))), 
                 intercept = lapply(1:ncolY, function(q)  rep(0, length = maxdim)), 
                 fitted    = lapply(1:ncolY, function(q)  matrix(0, nrow = nr, ncol = maxdim, dimnames = list(rownames(X), dimlab))), 
                 crit.reg  = rep(0, length = maxdim))
  names(res$XYcoef) <- names(res$intercept) <- names(res$fitted) <- colnames(Y)
  
  tabX <- X
  tabY <- as.data.frame(Y)
  lw   <-  ktabX$lw
  X.cw <- ktabX$cw
  blo  <- ktabX$blo
  rank <- maxdim
  eig  <- rep(0, maxdim)
  TL   <- ktabX$TL
  TC   <- ktabX$TC
  Yc1  <- matrix(0, nrow = ncolY, ncol = maxdim, dimnames = list(colnames(dudiY$tab), dimlab))
  l1   <- lY <- matrix(0, nrow = nr, ncol = maxdim, dimnames = list(row.names(dudiY$tab), dimlab))
  cov2 <- Ak <- matrix(0, nrow = nblo, ncol = maxdim, dimnames = list(names(ktabX$blo), dimlab))
  Tfa  <- lapply(1:nblo, function(k)  matrix(0, nrow = ncol(Xk[[k]]), ncol = maxdim, dimnames = list(colnames(Xk[[k]]), dimlab)))
  Tli  <- Tl1 <- rep(list(matrix(0, nrow = nr, ncol = maxdim, dimnames = list(row.names(dudiY$tab), dimlab))), nblo)
  faX  <- W   <- matrix(0, nrow = ncolX, ncol = maxdim, dimnames = list(col.names(ktabX), dimlab))
  
  
  
  ##-----------------------------------------------------------------------
  ##     Compute components and loadings by an iterative algorithm
  ##-----------------------------------------------------------------------
  
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  f1 <- function(x) lm.wfit(x = x, y = Y, w = lw)$fitted.values
  
  
  exp.mat <- function(MAT, EXP, tol=NULL){
    MAT    <- as.matrix(MAT)
    matdim <- dim(MAT)
    if(is.null(tol)){tol=min(1e-7, .Machine$double.eps*max(matdim)*max(MAT))}
    if(matdim[1]>=matdim[2]){ 
      svd1   <- svd(MAT)
      keep   <- which(svd1$d > tol)
      resmat <- t(svd1$u[,keep]%*%diag(svd1$d[keep]^EXP, nrow=length(keep))%*%t(svd1$v[,keep]))
    }
    if(matdim[1]<matdim[2]){ 
      svd1   <- svd(t(MAT))
      keep   <- which(svd1$d > tol)
      resmat <- svd1$u[,keep]%*%diag(svd1$d[keep]^EXP, nrow=length(keep))%*%t(svd1$v[,keep])
    }
    return(resmat)
  }
  
  
  for(h in 1 : maxdim) {
    
    ## Compute the matrix M for the eigenanalysis
    Gammak <- lapply(1:nblo, function(k) (1-gamma) * (t(Xk[[k]]) %*% Xk[[k]]) + gamma * diag(ktabX$blo[k]))
    Projk  <- lapply(1:nblo, function(k) Xk[[k]] %*% ginv(Gammak[[k]]) %*% t(Xk[[k]]))
    M      <- lapply(1:nblo, function(k) t(Y) %*% Projk[[k]] %*% Y)
    M      <- Reduce("+", M)
    
    ## Compute the loadings V and the components U (Y dataset)
    eig.M <- eigen(M)
    
    if (eig.M$values[1] < sqrt(.Machine$double.eps)) {
      rank <- h-1 ## update the rank
      break
    }
    
    eig[h]   <- eig.M$values[1]    
    Yc1[, h] <- eig.M$vectors[, 1, drop = FALSE]
    lY[, h]  <- Y %*% Yc1[, h]
    
    
    ## Compute the loadings Wk and the components Tk (Xk datasets)
    covutcarre <- 0
    covutk     <- rep(0, nblo)
    for (k in 1 : nblo) {
      Tfa[[k]][, h] <- (ginv(Gammak[[k]]) %*% t(Xk[[k]]) %*% lY[, h]) / sqrt(sum((exp.mat(Gammak[[k]], -0.5) %*% t(Xk[[k]]) %*% lY[, h])^2))
      Tl1[[k]][, h] <- Xk[[k]] %*% Tfa[[k]][, h]
      covutk[k]     <- crossprod(lY[, h] * lw, Tl1[[k]][, h])
      cov2[k, h]    <- covutk[k]^2  
      covutcarre    <- covutcarre + cov2[k, h]
    }
    
    for(k in 1 : nblo) {
      Ak[k, h]     <- covutk[k] / sqrt(sum(cov2[,h]))
      res$lX[, h]  <- res$lX[, h] + Ak[k, h] * Tl1[[k]][, h]
    }
    l1[, h] <- res$lX[, h] / sqrt(t(res$lX[, h])%*%res$lX[, h])
    W[, h]  <- tcrossprod(ginv(crossprod(X)), X) %*% res$lX[, h]
    
    ## Deflation of the Xk datasets on the global components T
    Xk <- lapply(Xk, function(y) lm.wfit(x = as.matrix(res$lX[, h]), y = y, w = lw)$residuals)
    X  <- as.matrix(cbind.data.frame(Xk))
  }
  
  
  ##-----------------------------------------------------------------------
  ##     Compute regressions coefficients
  ##-----------------------------------------------------------------------
  
  ## Use of the original (and not the deflated) datasets X and Y	
  X <- as.matrix(tabX)
  Y <- as.matrix(tabY)
  
  ## Computing the regression coefficients of X onto the global components T (Wstar)
  faX[, 1] <- W[, 1, drop = FALSE]
  A        <- diag(ncolX)
  if (maxdim >=2){
    for(h in 2:maxdim){
      a            <- t(l1[, h-1])%*%X / as.numeric((sqrt(t(res$lX[, h-1])%*%res$lX[, h-1])))
      A            <- A%*%(diag(ncolX) - W[, h-1]%*%a)
      faX[, h]     <- A%*%W[, h]
      X            <- X - l1[, h-1]%*%t(l1[, h-1])%*%X
    }
  }       
  
  ##  Compute the (eventually reduced or weighted) regression coefficients of X onto Y
  Yco           <- t(Y) %*% diag(lw) %*% res$lX
  norm.li       <- diag(crossprod(res$lX * sqrt(lw)))
  if (H==1){
    res$XYcoef <- lapply(1:ncolY, function(x) as.matrix(apply(sweep(faX, 2 , Yco[x,] / norm.li, "*"), 1, cumsum)))
  } else {
    res$XYcoef        <- lapply(1:ncolY, function(x) t(apply(sweep(faX, 2 , Yco[x,] / norm.li, "*"), 1, cumsum)))
  }
  names(res$XYcoef)    <- colnames(dudiY$tab)
  
  ## Correct the regression coefficients in case of block weghting (option == uniform)
  if (option[1] == "uniform"){
    In.X       <- unlist(sapply(1:nblo, function(k) rep(In.Xk[[k]], times=ktabX$blo[k])))
    res$XYcoef <- lapply(1:ncolY, function(x) res$XYcoef[[x]] * matrix(rep(In.Y, each=ncolX * maxdim), ncol = maxdim) / t(matrix(rep(In.X, each=maxdim), nrow=maxdim)))
  }
  
  
  ## Correct the regression coefficients in case of scaling (scale = TRUE)
  if (scale == TRUE){
    res$XYcoef <- lapply(1:ncolY, function(x) res$XYcoef[[x]] * matrix(rep(sdY[x], each=ncolX * maxdim), ncol = maxdim) / t(matrix(rep(sdX, each=maxdim), nrow=maxdim)))
  }
  
  ## Compute the intercept
  res$intercept        <- lapply(1:ncolY, function(x)  (meanY[x] - meanX %*% res$XYcoef[[x]]))
  names(res$intercept) <- colnames(dudiY$tab)
  
  
  ## Compute the regression error (/ raw data)
  rawX              <- as.matrix(cbind.data.frame(lapply(unclass(ktabX)[1 : nblo], scale, center = FALSE, scale = FALSE)))
  rawY              <- as.matrix(dudiY$tab)    
  res$fitted        <- lapply(1:ncolY, function(x) (matrix(rep(res$intercept[[x]], each=nr), nrow=nr) + rawX %*% res$XYcoef[[x]]))
  names(res$fitted) <- colnames(dudiY$tab)
  residual          <- lapply(1:ncolY, function(x) replicate(maxdim, rawY[, x]) - res$fitted[[x]])
  sum.residual.sq   <- lapply(1:ncolY, function(x) colSums(residual[[x]]^2))
  res$crit.reg      <- colSums(matrix(unlist(sum.residual.sq), nrow=ncolY, byrow = TRUE))/ncolY
  
  
  res$call      <- match.call()
  return(res)
}
