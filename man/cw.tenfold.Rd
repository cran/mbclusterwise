\name{cw.tenfold}


\alias{cw.tenfold}


\title{F-Fold cross-validation for clusterwise multiblock analyses}


\description{Function to perform a F-fold cross-validation applied to clusterwise multiblock analyses. This function is usually applied to various numbers of clusters and of dimensions to select their optimal values.}


\usage{cw.tenfold(Y, X, blo, option = c("none", "uniform"), G, H, FOLD = 10, INIT = 20, 
   method = c("mbpls", "mbpcaiv", "mbregular"), Gamma = NULL,
   parallel.level = c("high", "low"))}


\arguments{
  \item{Y}{a matrix or data frame containing the dependent variable(s)}
  \item{X}{a matrix or data frame containing the explanatory variables}
  \item{blo}{vector of the numbers of variables in each explanatory dataset}
  \item{option}{an option for the block weighting (by default, the first option is chosen): \cr
  	\option{none}{  the block weight is equal to the block inertia} \cr
  	\option{uniform}{  the block weight is equal to \eqn{1/K} for \eqn{(X_1, \dots, X_K)} and to 1 for \eqn{X} and \eqn{Y}}}
  \item{G}{an integer giving the number of clusters}
  \item{H}{an integer giving the number of dimensions of the component-based model}
  \item{FOLD}{an integer giving the number of folds of the F-Fold cross-validation procedure comprised between 2 and 10 (10 by default)}
  \item{INIT}{an integer giving the number of initializations required for the clusterwise analysis (20 by default)}
  \item{method}{an option for the multiblock method to be applied (by default, the first option is chosen):\cr
  	\option{mbpls}{  multiblock Partial Least Squares is applied}\cr
  	\option{mbpcaiv}{  multiblock Redundancy Analysis is applied}\cr
  	\option{mbregular}{  multiblock regularized regression is applied}}
  \item{Gamma}{a numeric value of the regularization parameter for the multiblock regularized regression comprised between 0 and 1 (NULL by default). The value (\code{Gamma=0}) leads to multiblock Redundancy Analysis and (\code{Gamma=1}) to multiblock PLS}
  \item{parallel.level}{Level of parallel computing, i.e. initializations are carried out simultaneously (high by default) \cr
 	\option{high}{  includes all the processing units of your computer}\cr
  	\option{low}{  includes only two processing units of your computer}}
}


\value{A list containing the following components is returned:
  \item{call}{the matching call}
  \item{sqrmse.cal}{the squared Root Mean Squared Error from the F calibration datasets}
  \item{sqrmse.val}{the squared Root Mean Squared Error from the F prediction datasets}
}


\references{Bougeard, S., Abdi, H., Saporta, G., Niang, N., Submitted, Clusterwise analysis for multiblock component methods.}


\author{Stephanie Bougeard (\email{stephanie.bougeard@anses.fr})}


\seealso{\code{\link{cw.multiblock}}, \code{\link{cw.predict}}}



\examples{
  data(simdata.red) 
  Data.X <- simdata.red[c(1:8, 21:28), 1:10]
  Data.Y <- simdata.red[c(1:8, 21:28), 11:13]
  res1   <- list()
  res2   <- list()

  ## Note that the options (INIT=2) and (parallel.level = "low") are chosen to quickly
  ## illustrate the function. 
  ## For real data, instead choose (INIT=20) to avoid local optima and (parallel.level = "high")
  ## to improve the computing speed.

  for (H in c(1:2)){
    print(paste("H=", H, sep=""))
    res1[[H]] <- cw.tenfold(Y = Data.Y, X = Data.X, blo = c(5, 5), option = "none", G = 1, H,  
      FOLD = 2, INIT = 2, method = "mbpls", Gamma = NULL, parallel.level = "low")
    res2[[H]] <- cw.tenfold(Y = Data.Y, X = Data.X, blo = c(5, 5), option = "none", G = 2, H,
       FOLD = 2, INIT = 2, method = "mbpls", Gamma = NULL, parallel.level = "low")
  }
  res1.cal <- unlist(lapply(1:2, function(x) mean(sqrt(res1[[x]]$sqrmse.cal), na.rm=TRUE))) 
  res1.val <- unlist(lapply(1:2, function(x) mean(sqrt(res1[[x]]$sqrmse.val), na.rm=TRUE)))
  res2.cal <- unlist(lapply(1:2, function(x) mean(sqrt(res2[[x]]$sqrmse.cal), na.rm=TRUE))) 
  res2.val <- unlist(lapply(1:2, function(x) mean(sqrt(res2[[x]]$sqrmse.val), na.rm=TRUE)))
  
  rmse.cal <- rbind(res1.cal, res2.cal)
  rmse.val <- rbind(res1.val, res2.val)
  rownames(rmse.cal) <- rownames(rmse.val) <- paste("G", 1:2, sep = "=")
  colnames(rmse.cal) <- colnames(rmse.val) <- paste("H", 1:2, sep = "=")

  par(mfrow=c(1,2))
  matplot(t(rmse.cal), type = "o", ylab = "RMSE of calibration", xlab = "Model dimension (H)", 
        main = "Calibration", col = c("steelblue", "darkorange"), pch = c(0, 5), lwd = c(3, 3))
        legend("center", inset = .05, legend = rownames(rmse.cal), pch = c(0, 5), lwd = c(3, 3),
        col = c("steelblue", "darkorange"), horiz = TRUE, title = "Cluster number (G)")
  matplot(t(rmse.val), type = "o", ylab = "RMSE of prediction", xlab = "Model dimension (H)", 
        main = "Prediction", col = c("steelblue", "darkorange"), pch = c(0, 5), lwd = c(3, 3))
        legend("center", inset = .05, legend = rownames(rmse.val), pch = c(0, 5), lwd = c(3, 3), 
        col = c("steelblue", "darkorange"), horiz = TRUE, title = "Cluster number (G)")
}


\keyword{multivariate}
\keyword{cluster}

