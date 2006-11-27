\name{svcm}
\alias{svcm}
\concept{sequential smoothing}
\concept{multidimensional smoothing}
\concept{tensor product}
\title{Fitting space-varying coefficient models}
\description{'svcm' is used to fit a 2d or 3d space-varying coefficient
  model or to merely smooth the input data using penalized
  B-splines. The smoothing works either sequentially or multidimensionally
  involving tensor products. So far, only space-invariant regressors are
  allowed. Data must be on a regular grid without missings.}  
\usage{
svcm(Y, X, vsize = c(1, 1, 1), knots = c(10, 10, 10),
     deg = c(1, 1, 1), opd = c(1, 1, 1), search = TRUE,
     lambda.init = rep(0.001, 3), method = "grid", type = "SEQ", ...)
%svcm(Y, X, vsize = c(1, 1, 1), knots = c(10, 10, 10), deg = c(1, 1, 1),
%     opd = c(1, 1, 1), search = TRUE, lambda.init = rep(0.001, 3), method =
%     "grid", type = "TP", taucs = FALSE, ...) 
}
\arguments{
  \item{Y}{array of observational data. Last dimension must correspond to
    the number of rows of \code{X}.}
  \item{X}{(r x p)-design matrix}
  \item{vsize}{numeric vector of the voxel size} 
  \item{knots}{vector of the numbers of inner knots in x-, y- (and z-)
    direction} 
  \item{deg}{vector of degrees of the basis functions in x-, y- (and z-)
    direction}
  \item{opd}{vector of the order of the difference penalties in x-, y-
    (and z-) direction}
  \item{search}{logical. If TRUE, the smoothing parameter will be
    optimized using \code{method} and GCV. If FALSE, \code{lambda.init}
    specifies the fixed smoothing parameter.}
  \item{lambda.init}{compulsory; initial value of global or dimension-specific
    smoothing parameter. See \bold{Details}.}
  \item{method}{optimization method to be used. See \bold{Details}.}
  \item{type}{character. \code{"SEQ"} (sequential) or \code{"TP"}
    (tensor product).} 
%  \item{taucs}{logical. If set TRUE, Taucs routines from Sivan
%    Toledo are employed. Use only if Taucs is installed on your
%    machine.}
  \item{\dots}{parameters to be passed to the optimization algorithm
    specified by \code{method}.}
}
\details{
  The purpose of \code{lambda.init} is three-fold: First, the length
  determines the use of either global or dimension-specific
  penalties. Second, it serves as fixed smoothing parameter if
  search is deactivated. Third, it is used as initial value from the
  \code{optim} algorithm which runs in case of a multidimensional tuning
  parameter when no grid search is desired.

  Unless method equals \code{"grid"}, \code{\link{optimize}} is called
  in the case of a global tuning parameter requiring a specified
  \code{interval} to be passed. While \code{\link{optimize}}
  does not take a starting value explicitly, a \code{startvalue} can be
  passed to \code{\link{cleversearch}}, e.g. \code{startvalue =
    lambda.init}.

  In the case of a dimension-specific tuning parameter, method
  \code{"grid"} evokes a full or greedy grid search (see 
  \code{\link{cleversearch}}). Amongst others, simplex method
  \code{"Nelder-Mead"} or quasi-Newton \code{"L-BFGS-B"} with positivity
  constraint for the smoothing parameter are conceivable, too. For further
  specification see \code{\link{optim}}.

  For simple smoothing of \code{Y} set \code{X = matrix(1, 1, 1)} and
  ascertain that the last dimension of \code{Y} matches \code{dim(X)[1]}.
}
\value{
  A list with components:
  \item{fitted}{fitted values as array of the same dimension as
    \code{Y}} 
  \item{effects}{effects of dimension (n.x, n.y, p) resp. (n.x, n.y,
    n.z, p).}
  \item{coeff}{coefficients (amplitudes of the basis functions) of
    dimension (p, r.x, r.y) resp. (p, r.x, r.y, r.z) with r.x number of
    basis functions in x-direction.}
  \item{knots}{see \code{knots}.}
  \item{deg}{see \code{deg}.}
  \item{opd}{see \code{opd}.}
  \item{vsize}{see \code{vsize}.}
  \item{type}{character describing the SVCM variant used. See \code{type}.}
  \item{call}{the matched call.}
  \item{opt}{a list with components depending on \code{search}, i.e. on
    whether optimization was performed or not:
    \describe{
      \item{time}{calculation/optimization time}
      \item{par}{initial value \code{lambda.init} or the best parameter
	found.}
      \item{value}{GCV value corresponding to \code{par}.}
      \item{GCVtab}{matrix of the search process with values of lambda
	            and corresponding GCV value.} 
      \item{...}{further values returned by optim(), optimize()}
    }
  }
}
\section{Warnings}{
  This model assumes the regressors to be space-invariant.
  Data must be on a regular grid without missings.
}
\note{The observations in \code{Y} are assigned to the center of the
  respective grid unit sized \code{vsize}. Hence the basis functions are
  evaluated at these coordinates.
}
\section{Background}{
  In the general case of 2d mesh data, Dierckx (1982, 1993) demonstrates
  the equivalence of successive univariate smoothing with smoothing
  based on a full bivariate B-spline matrix. However, the equivalence
  does no longer hold if penalties are introduced. Dierckx proposes the
  so-called 'new smoothing spline' as approximation to the
  multidimensional penalized smoothing (\code{type = "TP"}). While
  Dierckx determines the penalty structure through the spline degree and
  the equidistance between adjacent knots, the present implementation
  (\code{type = "SEQ"}) uses penalties of simple differences. 

  The calculation of GCV involves an inversion which is achieved using
  the recursion formula for band matrices in Hutchinson/de Hoog
  (1985). My collegue Thomas Kneib not only recommended this paper but
  also provided us with the basic.  
}
\references{
  Dierckx P. (1982) A fast algorithm for smoothing data on a rectangular
  grid while using spline functions. \emph{SIAM Journal on Numerical
    Analysis} \bold{19(6)}, 1286-1304.  

  Dierckx P. (1993) \emph{Curve and surface fitting with splines.}
  Oxford: Monographs on Numerical Analysis, Oxford University Press.

  Heim S., Fahrmeir L., Eilers P. H. C. and Marx B. D.
  (2006) Space-Varying Coefficient Models for Brain
  Imaging. \emph{Ludwig-Maximilians-University, SFB 386}, Discussion
  Paper \bold{455}.

  Hutchinson M. F. and de Hoog F. R. (1985) Smoothing noisy data with
  spline functions. \emph{Journal of Numerical Mathematics} \bold{47},
  99-106. 
  
%  Toledo S., Rotkin V. and Chen D. (2003) TAUCS: A library of sparse
%  linear solvers. Tel-Aviv University,
%  \url{http://www.tau.ac.il/~stoledo/taucs/} 
}
\seealso{
  \code{\link{optimize}} for one-dimensional minimization,\cr
  \code{\link{optim}} here explicitly method \code{"L-BGFS-B"},\cr
  \code{\link{cleversearch}} for clever or full grid search.
}
\examples{
## 2d DT-MRI data
data(brain2d)
X <- matrix(c(0.5, 0.5,   0,   0, 0.5, 0.5,
                0,   0, 0.5, 0.5, 0.5, 0.5,
              0.5, 0.5, 0.5, 0.5,   0,   0,
                0,   0,   0,   0,   1,  -1,
                1,  -1,   0,   0,   0,   0,
                0,   0,   1,  -1,   0,   0), nrow = 6)
##2d grid search for lambda; 60*50*6=18000 parameters (amplitudes) in total
M <- svcm(brain2d, X, knots=c(60, 50), deg=c(1, 1), vsize=c(1.875, 1.875),
          type="SEQ", lambda.init=rep(0.1, 2), search=TRUE,
          method="grid", lower=rep(-5, 2), upper=rep(0, 2), ngrid=10) 
str(M$opt)

## show results
zlim <- range(brain2d, M$fit)
old.par <- par(no.readonly=TRUE) 
par(pin=c(6, 1.2), mfrow=c(3, 1)) 
image(t(matrix(brain2d, nrow=dim(brain2d)[1])), axes=FALSE, zlim=zlim,
      col=grey.colors(256)) 
title("Observations: Six Diffusion Weighted Images")
image(t(matrix(M$fitted, nrow=dim(M$fit)[1])), axes=FALSE, zlim=zlim,
      col=grey.colors(256)) 
title("Fitted Values")
image(t(matrix(M$effects, nrow=dim(M$eff)[1])), axes=FALSE, 
      col=grey.colors(256))
title("Estimated Coefficient Surfaces: Six Elements of the Diffusion Tensor") 
par(old.par)


## 3d DT-MRI data; same regressors as in 2d; fixed global smoothing parameter
data(brain3d)
M3d <- svcm(brain3d, X, knots=c(5, 10, 5), deg=c(1, 1, 1), search=FALSE,
            vsize=c(1.875, 1.875, 4.0), type="TP", lambda.init=1)

## visualize results
zlim <- range(brain3d[,,,1], M3d$fit[,,,1])
old.par <- par(no.readonly = TRUE) 
par(pin=c(1.8, 5), mfrow=c(1, 3))
image(matrix(aperm(brain3d[,,,1], c(2,1,3)), nrow=dim(brain3d)[2]),
      axes=FALSE, zlim=zlim, col=grey.colors(256)) 
title("(a) Obs: 1st DWI")
image(matrix(aperm(M3d$fit[,,,1], c(2,1,3)), nrow=dim(brain3d)[2]),
      axes=FALSE, zlim=zlim, col=grey.colors(256)) 
title("(b) Fits of 1st DWI")
image(matrix(aperm(M3d$eff[,,,1], c(2,1,3)), nrow=dim(brain3d)[2]),
             axes=FALSE, col=grey.colors(256)) 
title("(c) Effects: 1st DT element")
title("Six axial slices of the 1st DWI-transform (a) and its fit (b);
      \n\n(c) corresponds to the first diffusion tensor component.",
       outer=TRUE, line=-5) 
par(old.par)
}
\keyword{smooth}   %Curve (and Surface) Smoothing
\keyword{spatial}  %Spatial Statistics
