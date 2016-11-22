\name{mbclusterwise-package}


\alias{mbclusterwise-package}


\alias{mbclusterwise}


\docType{package}


\title{\packageTitle{mbclusterwise}}


\description{\packageDescription{mbclusterwise}}


\details{
The DESCRIPTION file:
\packageDESCRIPTION{mbclusterwise}
\packageIndices{mbclusterwise}
}


\author{
\packageAuthor{mbclusterwise} \cr
Maintainer: \packageMaintainer{mbclusterwise}
}

\references{Bougeard, S., Abdi, H., Saporta, G., Niang, N., Submitted, Clusterwise analysis for multiblock component methods.}


\keyword{multivariate}
\keyword{cluster}


\seealso{\code{\link{ade4}}}


\examples{
  data(simdata.red) 
  Data.X <- simdata.red[c(1:10, 21:30), 1:10]
  Data.Y <- simdata.red[c(1:10, 21:30), 11:13]
  ## Note that the options (INIT=2) and (parallel.level = "low") are chosen to quickly
  ## illustrate the function. 
  ## For real data, instead choose (INIT=20) to avoid local optima and (parallel.level = "high")
  ## to improve the computing speed. 
  res.cw <- cw.multiblock(Y = Data.Y, X = Data.X, blo = c(5, 5), option = "none", G = 2, H = 1, 
            INIT = 2, method = "mbpls", Gamma = NULL, parallel.level = "low")
}
