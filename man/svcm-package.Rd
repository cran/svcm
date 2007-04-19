\name{svcm-package}
\alias{svcm-package}
\docType{package}
\title{
  2d and 3d Space-Varying Coefficient Models
}
\description{
  2d and 3d space-varying coefficient models are fitted to regular
  grid data using either a full B-spline tensor product approach or a
  sequential approximation. The latter one is computationally more
  efficient. Resolution increment is enabled.
}
\details{
  \tabular{ll}{
    Package: \tab svcm\cr
    Type: \tab Package\cr
    Version: \tab 0.1.1\cr
    Date: \tab 2007-04-19\cr
    License: \tab GPL version 2 or newer\cr
    Depends: \tab R (>= 2.4.0), Matrix, splines\cr
    LazyLoad: \tab yes\cr
    LazyData: \tab yes\cr
    URL: \tab http://www.statistik.lmu.de/~heim\cr
  }

  Originally, VCMs have been suggested by Hastie and Tibshirani (1993)
  for regressions with coefficients varying smoothly over a
  one-dimensional continuous variable such as time-varying effects.
  This package provides extensions to two- and three-dimensional
  space-varying coefficients surfaces for regularly gridded data without
  missings. Such a SVCM takes into account spatial correlation. The use of
  spline-basis functions serves to model the spatial coefficient
  field. As a consequence, estimates are accessible at any arbitrary
  position, not only on the original grid of voxels. Resolution can be
  easily increased and moreover penalized for possible initial
  anisotropy of the voxel size.  
  
  Two techniques are implemented. The multidimensional smoothing
  approach takes advantage of the sparsity of the spatial arrays
  involved. The second sequential one basically adapts the 'new
  smoothing spline' in Dierckx (1982), thus reducing the 3d (or
  higher-dimensional) problem to a sequence of one-dimensional
  smoothers. 

  The 2d and 3d examples have been chosen from the field of human brain
  imaging.
}
\author{
  Susanne Heim, with support from Paul Eilers, Thomas Kneib, and Michael Kobl

  Maintainer: Susanne Heim <susanneheim@gmx.net>
}
\references{
  Dierckx P. (1982) A fast algorithm for smoothing data on a rectangular
  grid while using spline functions. \emph{SIAM Journal on Numerical
    Analysis} \bold{19(6)}, 1286-1304.  

  Hastie T. and Tibshirani R. (1993) Varying-Coefficient Models (with
  discussion). \emph{Journal of the Royal Statistical Society} B
  \bold{55}, 757-796.

  Heim S., Fahrmeir L., Eilers P. H. C. and Marx B. D. (2006)
  Space-Varying Coefficient Models for Brain
  Imaging. \emph{Ludwig-Maximilians-University, SFB 386}, Discussion 
  Paper \bold{455}. 
}
\keyword{package}
