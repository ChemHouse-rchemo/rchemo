\name{consensuspca}
\alias{consensuspca}
\alias{transform.Consensuspca}
\alias{summary.Consensuspca}

\encoding{UTF-8}

\title{consensus PCA}

\description{

Algorithms fitting a consensus PCA of a list of matrices \eqn{Xlist}. 

A chosen PCA algorithm is applied on the X-matrix obtained after variable scaling, blockscaling, block concatenation.

\bold{Auxiliary functions}

\code{transform} Calculates the principal components for any new matrix \eqn{X} from the model.

\code{summary} returns summary information for the model.

}

\usage{

consensuspca(Xlist, blockscaling = TRUE, weights = NULL, nlv, 
Xscaling = c("none", "pareto", "sd")[1], 
algo = c("svd","eigen","eigenk","nipals","nipalsna","sph")[1], 
gs = TRUE, tol = .Machine$double.eps^0.5, maxit = 200)

\method{transform}{Consensuspca}(object, X, ..., nlv = NULL)  

\method{summary}{Consensuspca}(object, X, ...)  

}

\arguments{

\item{Xlist}{For the main function: list of training X-data (\eqn{n}rows).}

\item{X}{For the auxiliary functions: list of new X-data, with the same variables than the training X-data.}

\item{blockscaling}{logical. If TRUE, the scaling factor (computed on the training) is the "norm" of the block, i.e. the square root of the sum of the variances of each column of the block.}

\item{weights}{excepted for "nipalsna". Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{nlv}{For the main functions: The number(s) of PCs to calculate. --- For the auxiliary functions: The number(s) of PCs to consider.}

\item{Xscaling}{vector (of length Xlist) of variable scaling for each datablock, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{algo}{algorithm for the PCA among : "svd" (SVD factorization of \eqn{D^(1/2) * X}, using function \code{\link{svd}}), "eigen" (Eigen factorization of \eqn{X' * D * X}, using function \code{\link{eigen}}), "eigenk" (Eigen factorization of \eqn{D^(1/2) * X * X' D^(1/2)}, using function \code{\link{eigen}}), "nipals" (Eigen factorization of \eqn{X' * D * X} using NIPALS), "nipalsna" (Eigen factorization of \eqn{X' * D * X} using NIPALS allowing missing data in \eqn{X}), "sph" (Robust spherical PCA)}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

\bold{Specific for the NIPALS algorithm}

\item{gs}{Logical indicating if a Gram-Schmidt orthogonalization is implemented or not (default to \code{TRUE}).}

\item{tol}{Tolerance for testing convergence of the NIPALS iterations for each PC.}

\item{maxit}{Maximum number of NIPALS iterations for each PC.}


}

\value{

A list of outputs, such as

\item{T}{The X-score matrix (\eqn{n, nlv}).}

\item{Tk}{The partial X-score matrix from each Xlist matrix (\eqn{n, nlv}).}

\item{P}{The X-loadings matrix (\eqn{p, nlv}).}

\item{sv}{The singular values (\eqn{min(n, p), 1}) except for NIPALS = (\eqn{nlv, 1}).}

\item{eig}{The eigenvalues (\code{= sv^2}) (\eqn{min(n, p), 1}) except for NIPALS = (\eqn{nlv, 1}).}

\item{xmeans}{The list of centering vectors of \eqn{Xlist}.}

\item{xscales}{The list of \eqn{Xlist} variable standard deviations.}

\item{weights}{Weights applied to the training observations.}

\item{blockscaling}{block scaling.}

\item{Xnorms}{"norm" of each block, i.e. the square root of the sum of the variances of each column of each block, computed on the training, and used as scaling factor}.

\item{niter}{Numbers of iterations of the NIPALS.}

\item{conv}{Logical indicating if the NIPALS converged before reaching the maximal number of iterations.}

For \code{transform.Consensuspca}: X-scores matrix for new Xlist-data.

For \code{summary.Consensuspca}:

\item{explvarx}{matrix of explained variances.}

\item{contr.block}{block contributions.}

\item{contr.ind}{observation contributions.}

\item{contr.var}{variable contributions.}

\item{coord.var}{variable coordinates.}

\item{cor.circle}{variable coordinates on the correlation circle.}

}

\references{

Andrecut, M., 2009. Parallel GPU Implementation of Iterative PCA Algorithms. Journal of Computational Biology 16, 1593-1599. https://doi.org/10.1089/cmb.2008.0221

Gabriel, R. K., 2002. Le biplot - Outil d\'exploration de données multidimensionnelles. Journal de la Société Française de la Statistique, 143, 5-55.

Lingen, F.J., 2000. Efficient Gram-Schmidt orthonormalisation on parallel computers. Communications in Numerical Methods in Engineering 16, 57-66. https://doi.org/10.1002/(SICI)1099-0887(200001)16:1<57::AID-CNM320>3.0.CO;2-I

Tenenhaus, M., 1998. La régression PLS: théorie et pratique. Editions Technip, Paris, France.

Wright, K., 2018. Package nipals: Principal Components Analysis using NIPALS with Gram-Schmidt Orthogonalization. https://cran.r-project.org/

Wu, W., Massart, D.L., de Jong, S., 1997. The kernel PCA algorithms for wide data. Part I: Theory and algorithms. Chemometrics and Intelligent Laboratory Systems 36, 165-172. https://doi.org/10.1016/S0169-7439(97)00010-5


For Spherical PCA: 

Daszykowski, M., Kaczmarek, K., Vander Heyden, Y., Walczak, B., 2007. 
Robust statistics in data analysis - A review. Chemometrics and Intelligent 
Laboratory Systems 85, 203-219. https://doi.org/10.1016/j.chemolab.2006.06.016

Locantore N., Marron J.S., Simpson D.G., Tripoli N., Zhang J.T., Cohen K.L.
Robust principal component analysis for functional data, Test 8 (1999) 1-7

Maronna, R., 2005. Principal components and orthogonal regression based on 
robust scales, Technometrics, 47:3, 264-273, DOI: 10.1198/004017005000000166

}


\examples{

n <- 10 ; p <- 10
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- rnorm(n)

m <- 2
Xtest <- matrix(rnorm(m * p), ncol = p)

colnames(Xtrain) <- colnames(Xtest) <- paste("v", 1:p, sep = "")

Xtrain
Xtest

blocks <- list(1:2, 4, 6:8)
X1 <- mblocks(Xtrain, blocks = blocks)
X2 <- mblocks(Xtest, blocks = blocks)

nlv <- 3
fm <- consensuspca(Xlist = X1, Xscaling = c("sd","none","none"), 
blockscaling = TRUE, weights = NULL, nlv = nlv)

summary(fm, X1)
transform(fm, X2)

}

\keyword{datagen}