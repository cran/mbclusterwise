\name{cw.multiblock}


\alias{cw.multiblock}


\title{Clusterwise multiblock analyses}


\description{Function to perform a clusterwise multiblock analyses (clusterwise multiblock Partial Least Squares, clusterwise multiblock Redundancy Analysis or clusterwise regularized multiblock regression) of several explanatory blocks \eqn{(X_1, \dots, X_K)} to explain a dependent dataset \eqn{Y}.}


\usage{cw.multiblock(Y, X, blo, option = c("none", "uniform"), G, H, INIT = 20, 
    method = c("mbpls", "mbpcaiv", "mbregular"), Gamma = NULL, 
    parallel.level = c("high", "low"))}



\arguments{
  \item{Y}{a matrix or data frame containing the dependent variable(s)}
  \item{X}{a matrix or data frame containing the explanatory variables}
  \item{blo}{a vector of the numbers of variables in each explanatory dataset}
  \item{option}{an option for the block weighting (by default, the first option is chosen): \cr
  	\option{none}{  the block weight is equal to the block inertia} \cr
  	\option{uniform}{  the block weight is equal to \eqn{1/K} for \eqn{(X_1, \dots, X_K)} and to 1 for \eqn{X} and \eqn{Y}}}
  \item{G}{an integer giving the expected number of clusters}
  \item{H}{an integer giving the expected number of dimensions of the component-based model}
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
  \item{error}{a vector containing the value of the criterion to be minimized (overall prediction error) ; this error is performed on the centered and scaled data}
  \item{beta.cr}{a list of array that contain the intercept and the regression coefficients associated with the centered and scaled data for each of the \eqn{G} clusters}
  \item{beta.raw}{a list of array that contain the intercept and the regression coefficients associated with the raw data for each of the \eqn{G} clusters}
  \item{hopt}{the real number of dimensions of the component-based model (\code{hopt} is sometimes lower than the expected \code{H})}
  \item{Ypred.cr}{a list of matrices that contain the predicted dependent values associated with the centered and scaled data for each of the \eqn{G} clusters}
  \item{Ypred.raw}{a list of matrices that contain the predicted dependent values associated with the raw data for each of the \eqn{G} clusters}
  \item{cluster}{a vector containing the observation assignation to the \eqn{G} expected clusters (when \eqn{G}>1 only)}
}


\references{Bougeard, S., Abdi, H., Saporta, G., Niang, N., Submitted, Clusterwise analysis for multiblock component methods.}


\author{Stephanie Bougeard (\email{stephanie.bougeard@anses.fr})}


\seealso{\code{\link{cw.tenfold}}, \code{\link{cw.predict}}}



\examples{
  data(simdata.red) 
  Data.X <- simdata.red[c(1:10, 21:30), 1:10]
  Data.Y <- simdata.red[c(1:10, 21:30), 11:13]
  ## Note that the options (INIT=2) and (parallel.level = "low") are chosen to quickly
  ## illustrate the function. 
  ## For real data, instead choose (INIT=20) to avoid local optima and (parallel.level = "high")
  ## to improve the computing speed.
  res.cw <- cw.multiblock(Y = Data.Y, X = Data.X, blo = c(5, 5), option = "none", G = 2, 
            H = 1, INIT = 2, method = "mbpcaiv", Gamma = NULL, parallel.level = "low")
}


\keyword{multivariate}
\keyword{cluster}


