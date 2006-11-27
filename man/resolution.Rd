\name{resolution}
\alias{resolution}
\title{Re-scaling resolution of SVCM predictors and effects}
\description{
  This routine serves post-hoc adjustment of the resolution of a
  space-varying coefficient model estimated by \code{\link{svcm}}.
}
\usage{
resolution(X, svcmlist, fac)
}
\arguments{
  \item{X}{(r x p)-array of covariates}
  \item{svcmlist}{return list of function \code{\link{svcm}}}
  \item{fac}{2d or 3d vector of scaling factors} 
}
\details{
  The basis function approach underlying \code{\link{svcm}} allows to
  rescale the original resolution by evaluating the basis functions at
  additional points. Assuming that the voxel center is most
  representative for the whole voxel, \code{fac}-times resolution of 1d
  data with \eqn{n} voxels sized \eqn{vsize} bases on coordinates
  \deqn{\left(i - \frac{1}{2}\right) \cdot \frac{vsize}{fac} , \quad i = 1,
    \ldots, n \cdot fac.}{(i - 0.5) * vsize/fac, i = 1, ..., n*fac.}  
  See also doc file resolution\_scheme.pdf.\cr
  The formula is applied into x-, y- and z-direction and results in a
  refined equidistant 2d resp. 3d grid. It also means that, in general,
  the resized arrays of predictors and effects do no longer contain the
  values at the former coordinates.

  Note that memory requirements can be enormous depending on object size
  and the intended resolution. 
}
\value{
  A list with components  
  \item{fitted}{fitted values at \code{fac}-times rescaled resolution.}
  \item{effects}{estimated effects at \code{fac}-times rescaled resolution.}
}
\examples{
##DTI data; regressors are given by the diffusion weigthing gradients
data(brain2d)
X <- matrix(c(0.5, 0.5,   0,   0, 0.5, 0.5,
                0,   0, 0.5, 0.5, 0.5, 0.5,
              0.5, 0.5, 0.5, 0.5,   0,   0,
                0,   0,   0,   0,   1,  -1,
                1,  -1,   0,   0,   0,   0,
                0,   0,   1,  -1,   0,   0), nrow = 6)
M <- svcm(brain2d, X, knots=c(60, 50), deg=c(1, 1), vsize=c(1.875,
          1.875), search=TRUE, type="SEQ", lambda.init=rep(0.1, 2),
          lower=rep(-5, 2), upper=rep(0, 2), ngrid=10)
M2 <- resolution(X, M, fac=c(2, 2))

##show data extract at original and double resolution
extract <- list(M$fit[21:40, 21:60, 1],
                M2$fit[(21*2):(40*2), (21*2):(60*2), 1], 
                M$eff[21:40, 21:60, 1],
                M2$eff[(21*2):(40*2), (21*2):(60*2), 1]) 
zlim1 <- range(extract[[1]], extract[[2]])
zlim2 <- range(extract[[3]], extract[[4]])
old.par <- par(no.readonly = TRUE) 
par(pin=c(3*1, 3*0.67), mfrow=c(2, 2)) 
image(t(extract[[1]]), axes=FALSE, zlim=zlim1, col=grey.colors(256))
title("Fitted Values")
image(t(extract[[2]]), axes=FALSE, zlim=zlim1, col=grey.colors(256))
title("Fitted Values at Double Resolution")
image(t(extract[[3]]), axes=FALSE, zlim=zlim2, col=grey.colors(256))
title("Estimated VC Surface (1st DT element)")
image(t(extract[[4]]), axes=FALSE, zlim=zlim2, col=grey.colors(256))
title("VC Surface at Double Resolution")
par(old.par)
}
\keyword{smooth}
\keyword{spatial}
