\name{covselrda}
\alias{covselrda}
\alias{covsellda}
\alias{covselqda}
\alias{predict.Covselrda}
\alias{predict.Covselprobda}
\encoding{latin1}

\title{CovSel-discriminant analysis}

\description{

Variable selection for high-dimensionnal data with the COVSEL method (Roger et al. 2011), followed by a linear regression model.

The training variable \eqn{y} (univariate class membership) is firstly transformed to a dummy table containing \eqn{nclas} columns, where \eqn{nclas} is the number of classes present in \eqn{y}. Each column is a dummy variable (0/1). Then, a variable selection, based on the COVSEL method, is implemented on the \eqn{X-}data and the dummy table, returning a set of X-variables that are used as dependent variables in a DA model.

- \code{covselrda}: A linear regression model predicts the Y-dummy table from the selected X-variables. For a given observation, the final prediction is the class corresponding to the dummy variable for which the prediction is the highest.

- \code{covsellda} and \code{covselqda}: Probabilistic LDA and QDA are run over the selected X-variables, respectively.


\bold{Auxiliary functions}

\code{predict} Calculates the predictions for any new set of variables contained in the selection.

}

\usage{

covselrda(X, y, nvar = NULL, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL)

covsellda(X, y, nvar = NULL, prior = c("unif", "prop"), Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL)
  
covselqda(X, y, nvar = NULL, prior = c("unif", "prop"), Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL)

predict.Covselrda(object, X, ..., nvar = NULL)

predict.Covselprobda(object, X, ..., nvar = NULL)

}

\arguments{

\item{X}{X-data (\eqn{n, p}).}

\item{y}{Training class membership (\eqn{n}). \bold{Note:} If \code{y} is a factor, it is replaced by a character vector.}

\item{nvar}{Number of variables to select in \code{X}. Can be a vector for the auxiliary functions}

\item{prior}{The prior probabilities of the classes. Possible values are "unif" (default; probabilities are set equal for all the classes) or "prop" (probabilities are set equal to the observed proportions of the classes in \code{y}).}

\item{Xscaling}{X variable scaling among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{Yscaling}{Y variable scaling among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{weights}{Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

}

\value{

For \code{covselrda}, \code{covselda}, \code{covselqda}:

\item{sel}{A dataframe where variable \code{sel} shows the column indexes of the variables selected in \eqn{X}.}

\item{fm}{List of linear regression or discriminant models, involving 1 to nvar selected explicative variables.}

\item{lev}{classes}

\item{ni}{number of observations in each class}

\item{weights}{The weights used for the row observations.}


For \code{predict.Covselrda}, \code{predict.Covselprobda}: 

\item{pred}{predicted class for each observation}

\item{posterior}{calculated probability of belonging to a class for each observation}

}

\references{

Roger, J.M., Palagos, B., Bertrand, D., Fernandez-Ahumada, E., 2011. CovSel: Variable selection for highly multivariate and multi-response calibration: Application to IR spectroscopy. Chem. Lab. Int. Syst. 106, 216-223. 

}


\examples{

n <- 6 ; p <- 4
X <- matrix(rnorm(n * p), ncol = p)
y <- c("A","A","B","B","C","C")

sel <- covselrda(X, y, nvar = 3)

predict(sel, X, nvar = c(2,3))

}

\keyword{datagen}