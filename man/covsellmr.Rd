\name{covsellmr}
\alias{covsellmr}
\alias{predict.Covsellmr}
\encoding{latin1}

\title{CovSel-linear regression model}

\description{

Variable selection for high-dimensionnal data with the COVSEL method (Roger et al. 2011), followed by a linear regression model.

\bold{Auxiliary functions}

\code{predict} Calculates the predictions from the regression model for any new set of variables contained in the selection.

}

\usage{

covsellmr(X, Y, nvar = NULL, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL)

predict.Covsellmr(object, X, ..., nvar = NULL)

}

\arguments{

\item{X}{X-data (\eqn{n, p}).}

\item{Y}{Y-data (\eqn{n, q}).}

\item{nvar}{Number of variables to select in \code{X}.}

\item{Xscaling}{X variable scaling among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{Yscaling}{Y variable scaling among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{weights}{Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

}

\value{

\item{sel}{A dataframe where variable \code{sel} shows the column indexes of the variables selected in \eqn{X}.}

\item{fm}{List of linear regression models, involving 1 to nvar selected explicative variables.}

\item{weights}{The weights used for the row observations.}

}

\references{

Roger, J.M., Palagos, B., Bertrand, D., Fernandez-Ahumada, E., 2011. CovSel: Variable selection for highly multivariate and multi-response calibration: Application to IR spectroscopy. Chem. Lab. Int. Syst. 106, 216-223. 

}


\examples{

n <- 6 ; p <- 4
X <- matrix(rnorm(n * p), ncol = p)
Y <- matrix(rnorm(n * 2), ncol = 2)

sel <- covsellmr(X, Y, nvar = 3)

predict(sel, X, nvar = c(2,3))

}

\keyword{datagen}