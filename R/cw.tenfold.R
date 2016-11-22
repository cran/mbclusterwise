cw.tenfold <-
function(Y, X, blo, option = c("none", "uniform"), G, H, FOLD = 10, INIT = 20, method = c("mbpls", "mbpcaiv", "mbregular"), Gamma = NULL, parallel.level = c("high", "low")){
  
  

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
  if (FOLD < 2)
    stop("Minimum two folds are required")
  if (FOLD > 10)
    stop("Maximum ten folds are required")
  if (INIT < 2)
    stop("Minimum two initializations are required")
  if (method == "mbpls" | method == "mbpcaiv"){
    if (!is.null(Gamma)) stop("Non convenient value for parameter Gamma")
  } else {
    if (is.null(Gamma)) stop("Non convenient value for parameter Gamma")
  }
  option <- match.arg(option)
  method <- match.arg(method)
  
  
  
  # ---------------------------------------------------------------------------
  # 1. Required parameters
  # ---------------------------------------------------------------------------
  
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  Q <- dim(Y)[2]
  
  
  # ---------------------------------------------------------------------------
  # 2. Randomly shuffle the data by row 
  # ---------------------------------------------------------------------------
  
  set.seed((as.numeric(Sys.time()) - floor(as.numeric(Sys.time()))) * 1e8 -> seed)
  sample.obs <- sample(nrow(X))
  X.sample   <- X[sample.obs, ]
  Y.sample   <- Y[sample.obs, , drop = FALSE]
  folds      <- cut(seq(1, dim(X)[1]), breaks = FOLD, labels = FALSE)
  
  
  # ---------------------------------------------------------------------------
  # 3. Apply the cross-validation procedure
  # ---------------------------------------------------------------------------
  
  sqrmse.cal <- c()
  sqrmse.val <- c()
  
  for (i in 1 : FOLD){
    
    print(paste("Fold", i))
    
    # ---------------------------------------------------------------------------
    # 31. Create datasets of equal size to get train and test data
    # ---------------------------------------------------------------------------    
    
    test.row <- which(folds == i, arr.ind = TRUE)
    X.test   <- X.sample[test.row, ]
    X.train  <- X.sample[-test.row, ]
    Y.test   <- Y.sample[test.row, , drop = FALSE]
    Y.train  <- Y.sample[-test.row, , drop = FALSE]
    
    
    # ---------------------------------------------------------------------------
    # 32. Apply the clusterwise procedure (model & prediction)
    # ---------------------------------------------------------------------------   
    
    # Process the clusterwise regression with the train datasets
    rescw.train <- cw.multiblock(Y = Y.train, X = X.train, blo, option, G, H, INIT, method, Gamma, parallel.level) 

    # Process the X.test prediction from the train model
    rescw.pred  <- cw.predict(Xnew = X.test, res.cw = rescw.train)

    
    # ---------------------------------------------------------------------------
    # 33. Evaluate the clusterwise regression performance
    # ---------------------------------------------------------------------------   
    
    # Comparison of the observed and predicted Y.train
    sqrmse.cal[i] <- min(rescw.train$error)
    
    # Comparison of the observed and predicted Y.test (normalized according to Y.train)
    Y.test.c  <- sweep(Y.test, 2, colMeans(Y.train), FUN = "-")
    Y.test.cr <- sweep(Y.test.c, 2, apply(Y.train, 2, sd), FUN = "/")
    if (G == 1){
      sqrmse.val[i] <- sum((Y.test.cr - rescw.pred$Ypred.cr)^2) / Q
    } else {
      G.pred        <- as.numeric(levels(factor(rescw.pred$clusternew)))
      error.g       <- sapply(G.pred, function(g) sum((Y.test.cr[which(rescw.pred$clusternew == g), , drop = FALSE] - rescw.pred$Ypred.cr[[g]])^2) / Q)
      sqrmse.val[i] <- sum(error.g)
    }
    
  }
  
  
  # ---------------------------------------------------------------------------
  #  4. Result storage
  # ---------------------------------------------------------------------------
  
  res            <- list()
  res$sqrmse.cal <- sqrmse.cal
  res$sqrmse.val <- sqrmse.val
  return(res) 
}
